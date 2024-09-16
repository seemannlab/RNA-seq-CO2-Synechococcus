# Make a summary plot for the Cluster enrichment

library(tidyverse)
library(ggpubr)

library(treemapify)

################################################################################
# Load input

annot <-
  'data/C_annotation.tsv' |>
  read_tsv()

deg <-
  'analysis/D_stagewise-adjusted-DEGs.tsv' |>
  read_tsv()

vst <-
  'analysis/D_vst-expression.tsv' |>
  read_tsv()

meta <-
  'data/C_meta.tsv' |>
  read_tsv() |>
  mutate_at('CO2', ~ fct_reorder(as.character(.x), as.numeric(.x)))


################################################################################
# Z-scale VST expression with subsequent averaging

vst.mat <-
  vst |>
  select(- Geneid) |>
  as.matrix() |>
  magrittr::set_rownames(vst$Geneid)

z.mat <-
  vst.mat |>
  apply(1, scale) |>
  t() |> 
  magrittr::set_colnames(colnames(vst.mat))

z.expr <-
  z.mat |>
  as_tibble(rownames = 'Geneid') |>
  pivot_longer(- Geneid, names_to = 'lib', values_to = 'z.expression') |>
  left_join(meta, 'lib')

z.avg <-
  z.expr |>
  group_by(Geneid, CO2) |>
  summarize(avg.z = mean(z.expression)) |>
  ungroup() |>
  # add locus tag with taxid for stringapp (later use)
  left_join(annot, 'Geneid') |>
  mutate(string = ifelse(
    !is.na(old_locus_tag),
    paste0('32049.', old_locus_tag),
    NA_character_
  ))

write_tsv(z.avg, 'analysis/E_zvst.tsv')

################################################################################
# download all brite hierarchy data from KEGG

kegg.brite <-
  'https://rest.kegg.jp/link/syp/brite' |>
  read_tsv(col_names = c('brite', 'syp'))

kegg.brite.json <-
  kegg.brite |>
  pull(brite) |>
  unique() %>%
  paste0('https://rest.kegg.jp/get/', ., '/json') |>
  map(read_lines)

kegg.brite.names <-
  'https://rest.kegg.jp/list/br' |>
  read_tsv(col_names = c('br', 'name')) |>
  with(set_names(
    name,
    paste0('syp', br)
  ))

################################################################################
# Convert KEGG hierarchy into a nice table

# xs <- jsonlite::parse_json(kegg.brite.json[[1]] )

build_table <- function(xs, i = 1, pre = NULL) {
  if ('children' %in% names(xs)) {
    # can do recursion
    y <- paste0('level', i)
    pre2 <- bind_cols(pre, tibble(!!y := xs$name))
    xs$children |>
      map(build_table, i = i + 1, pre = pre2) |>
      bind_rows()
    
  } else {
    # return accumulated results
    y <- paste0('level', i)
    bind_cols(pre, tibble(!!y := xs$name))
  }
}

brite.tbl <-
  kegg.brite.json |>
  map(safely(jsonlite::parse_json)) |>
  map('result') |>
  map(build_table) |>
  bind_rows() |>
  mutate(row = 1:n())
  
brite.tbl.long <-
  brite.tbl |>
  pivot_longer(- row)

# tidy up with gene names separately in one column
brite.genes <-
  brite.tbl.long |>
  filter(str_detect(value, 'SYNPCC7002_')) |>
  select(row, gene = value) |>
  separate_rows(gene, sep = '\t') |>
  filter(str_detect(gene, 'SYNPCC7002_')) |>
  mutate(gene = str_remove(gene, ' .*$')) |>
  rename(old_locus_tag = gene)


# tidy up the levels
brite.info <-
  brite.tbl.long |>
  drop_na() |>
  filter(! str_detect(value, 'SYNPCC7002_')) |>
  # remove all eukaryotic rows
  anti_join(
    brite.tbl.long |>
      filter(value == 'Eukaryotic type'),
    'row'
  ) |>
  # remove notes on prokatyoric entries (but keep the rows)
  filter(value != 'Prokaryotic type') |>
  # fix leveling where needed (after removing some entries)
  group_by(row) |>
  mutate(name = paste0('level', 1:n())) |>
  spread(name, value, fill = '')

# combine into one large table
brite <-
  brite.genes |>
  inner_join(brite.info, 'row') |>
  mutate(level1 = kegg.brite.names[level1]) |>
  select(- row)

write_tsv(brite, 'analysis/E_brite-hierarchy.tsv')
  
################################################################################

z.avg |>
  transmute(
    old_locus_tag,
    CO2 = paste0(CO2, '%') |>
      fct_reorder(CO2 |> as.numeric()),
    avg.z
  ) |>
  inner_join(brite, 'old_locus_tag') |>
  filter(
    level1 == 'KEGG Orthology (KO)',
    level2 != '09190 Not Included in Pathway or Brite',
    level2 != '09180 Brite Hierarchies'
  ) |>
  mutate_at(vars(contains('level')), str_remove, ' \\[.*$') |>
  mutate_at(vars(contains('level')), str_remove, '^[0-9]+ ') -> brite.ko


brite.ko |>
  semi_join(
    deg |>
      filter(is.de) |>
      select(old_locus_tag) |>
      unique(),
    'old_locus_tag'
  ) |>
  mutate(
    CO2 = CO2 |>
      paste('CO2') |>
      fct_relevel(
        '0.04% CO2',
        '4% CO2',
        '30% CO2',
        '8% CO2'
      )
  ) |>
  ggplot(aes(
    area = 1, fill = avg.z,
    subgroup = level2,
    subgroup2 = level3,
    subgroup3 = level4
  )) +
  scale_fill_gradient2(
    low =  RColorBrewer::brewer.pal(3, 'RdBu')[3],
    high =  RColorBrewer::brewer.pal(3, 'RdBu')[1],
    name = 'z-scaled expression'
  ) +
  geom_treemap() +
  geom_treemap_subgroup_border(color = 'black', size = 10) +
  geom_treemap_subgroup2_border(colour = "black", size = 5) +
  geom_treemap_subgroup3_border(colour = "black", size = 1) +
  # geom_treemap_subgroup2_text(
  #   place = "centre", grow = T, alpha = 0.9,
  #   reflow = TRUE,
  #   colour = "yellow",
  #   min.size = 2
  # ) +
  geom_treemap_subgroup3_text(
    place = "centre", grow = T, alpha = 0.9,
    reflow = TRUE,
    colour = "black",
    min.size = 2
  ) +
  facet_wrap(~ CO2) +
  theme_pubr(18) +
  theme(legend.key.width = unit(3, 'cm'))

ggsave('analysis/E_treemap.jpeg', width = 16, height = 12, dpi = 500)

################################################################################
################################################################################
sessionInfo()
