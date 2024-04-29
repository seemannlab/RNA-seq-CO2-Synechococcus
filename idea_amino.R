#Purpose: Check general expression patterns in Amino Acid biosynthesis

library(tidyverse)

# All SYP7002 AA metabolic pathways
paths <-
  tribble(
    ~ kegg, ~ title,
    'syp00250', 'Alanine, aspartate and glutamate metabolism',
    'syp00260', 'Glycine, serine and threonine metabolism',
    'syp00270', 'Cysteine and methionine metabolism',
    'syp00280', 'Valine, leucine and isoleucine degradation',
    'syp00290', 'Valine, leucine and isoleucine biosynthesis',
    'syp00300', 'Lysine biosynthesis',
    'syp00310', 'Lysine degradation',
    'syp00220', 'Arginine biosynthesis',
    'syp00330', 'Arginine and proline metabolism',
    'syp00340', 'Histidine metabolism',
    'syp00350', 'Tyrosine metabolism',
    'syp00360', 'Phenylalanine metabolism',
    'syp00380', 'Tryptophan metabolism',
    'syp00400', 'Phenylalanine, tyrosine and tryptophan biosynthesis'
)

################################################################################


paths.gene <-
  'analysis/I_gene2pathway.tsv' |>
  read_tsv() |>
  semi_join(paths, c('Pathway' = 'kegg'))

paths.gene |>
  plyr::dlply('Title', \(x) x$old_locus_tag) |>
  venn::venn(zcolor = 'style', ilcs = 2, sncs = 2)

################################################################################

library(ggkegg)

################################################################################

annot <-
  'data/C_annotation.tsv' |>
  read_tsv()

vst <-
  'analysis/D_vst-expression.tsv' |>
  read_tsv()

meta <-
  'data/C_meta.tsv' |>
  read_tsv() |>
  mutate_at('CO2', ~ fct_reorder(as.character(.x), as.numeric(.x)))

################################################################################
# z-scale expression

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

################################################################################
# Test for CyanoCyc

deg <-
  'analysis/D_stagewise-adjusted-DEGs.tsv' |>
  read_tsv()

deg |>
  filter(is.de, abs(log2FoldChange) >= 1) |>
  select(Geneid) |>
  unique() -> de.mask

z.expr |>
  filter(str_detect(Geneid, 'gene')) |>
  group_by(Geneid, CO2) |>
  summarize(avg = mean(z.expression)) |>
  ungroup() |>
  mutate(CO2 = paste0(CO2, 'pct')) |>
  spread(CO2, avg) |>
  semi_join(de.mask, 'Geneid') |>
  mutate(Geneid = str_remove(Geneid, 'gene-')) |>
  select(Geneid, `0.04pct`, `4pct`, `8pct`, `30pct`) |>
  write_tsv('~/Downloads/cyanocyc-vst-scaled.tsv')

expression <-
  'analysis/D_normalized-counts.tsv' |>
  read_tsv()

expression |>
  filter(str_detect(Geneid, 'gene')) |>
  pivot_longer(- Geneid, names_to = 'lib') |>
  left_join(meta, 'lib') |>
  group_by(Geneid, CO2) |>
  summarize(avg = mean(value)) |>
  ungroup() |>
  mutate(CO2 = paste0(CO2, 'pct')) |>
  spread(CO2, avg) |>
  semi_join(de.mask, 'Geneid') |>
  mutate(Geneid = str_remove(Geneid, 'gene-')) |>
  select(Geneid, `0.04pct`, `4pct`, `8pct`, `30pct`) |>
  write_tsv('~/Downloads/cyanocyc-sf.tsv')
           
  

################################################################################
# Adaptation of function from `L_pathways` to also show metabolic pathways
# without explicit gene positions in figure
# Thus: Extract the average x/y position of the arrow ploygon's coords
# plot heatmap boxes there

my.extra.path <- function(x, legend.pos) {
  # x <- 'syp01230'
  # Load pathway
  path <-
    x |>
    pathway()
  
  # extract genes from pathway
  path.genes <-
    path |>
    activate(nodes) |>
    as_tibble() |>
    filter(type == 'gene')
  
  # join to expression
  path.expr <-
    path.genes |>
    select(name) |>
    unique() |>
    # match locus to Geneid
    mutate(locus = str_remove_all(name, 'syp:')) |>
    separate_rows(locus, sep = ' ') |>
    left_join(annot, c('locus' = 'old_locus_tag')) |>
    select(name = name.x, Geneid) |>
    # add z-scaled expression
    left_join(z.expr, 'Geneid', relationship = "many-to-many") |>
    # average per 'enzyme' NOT gene (could be multiple)
    group_by(name, CO2) |>
    summarize(avg = mean(z.expression, na.rm = TRUE)) |>
    ungroup() |>
    # left_join(path.genes, 'name', relationship = "many-to-many") |>
    # ensure all values fit the ±2 scale
    mutate(avg = pmin(pmax(avg, -2), 2))
  
  
  # Extract polygon of gene reaction
  path |>
    activate(nodes) |>
    as_tibble() |>
    filter(type == 'gene') |>
    select(orig.id, name, reaction, graphics_name, coords) |>
    #! symbol '|' should indicate a new line/polygon for a "branch"
    # regex fits ',|' between two number pairs that are split per ','
    mutate_at('coords', str_extract_all, '(?<=[,|]?)\\d+,\\d+') |>
    unnest(coords) |>
    separate(coords, c('x', 'y'), sep = ',') |>
    mutate_at(c('x', 'y'), as.integer) |>
    # find 'center' point
    group_by(orig.id, name, reaction, graphics_name) |>
    summarize(x = mean(x), y = - mean(y)) |>
    ungroup() -> gene.centers
    
  
  # base KEGG path viz
  p.path <-
    path |>
    ggraph(layout="manual", x=x, y=y) +
    # ensure the centers look right
    # geom_point(
    #   aes(x, y),
    #   size = 5, color = 'yellow', alpha = .5,
    #   data = gene.centers
    # ) +
    overlay_raw_map() +
    theme_void()
  # p.path
  
  # add expression for genes into the base visualization
  dims <- list(w = 20, h = 10)
  helper <- function(gene_row) {
    p <-
      gene_row |>
      as_tibble() |>
      select(name) |>
      left_join(path.expr, 'name', relationship = 'one-to-many') |>
      ggplot(aes(CO2, 1, fill = avg)) +
      geom_tile(color = 'black') +
      # coord_fixed() +
      coord_polar() +
      scale_fill_gradientn(colors = RColorBrewer::brewer.pal(5, 'RdBu') |> rev(),
                           # make sure all have the same scale
                           limits = c(-2, 2)) +
      theme_void() +
      theme(legend.position = 'hide')
    # helper annotation to add to main plot
    gene_row |>
      with(annotation_custom(
        ggplotGrob(p), 
        xmin = x - dims$w, xmax = x + dims$w,
        ymin = y - dims$h, ymax = y + dims$h
      ))
  }
  
  gene.plots <-
    gene.centers |>
    rowwise() |>
    do(i = helper(.)) |>
    pull(i)
  
  
  # create legend
  p.legend <-
    path.expr |>
    # arbitrarily choose first one for color bar
    filter(name == first(name)) |>
    ggplot(aes(CO2, 1, fill = avg, label = CO2)) +
    geom_tile() +
    coord_polar() +
    scale_fill_gradientn(colors = RColorBrewer::brewer.pal(5, 'RdBu') |> rev(),
                         name = 'Average z-scaled expression',
                         limits = c(-2, 2)) +
    guides(fill = guide_colorbar(title.position = 'top', barwidth = unit(5, 'cm'))) +
    ggpubr::theme_pubclean(8)  +
    theme(
      axis.title.y = element_blank(),
      axis.text.y =  element_blank(),
      axis.ticks = element_blank(),
      panel.grid = element_blank(),
      panel.grid.major.y = element_blank()
    ) +
    geom_tile(fill = I('white'), color = 'black') +
    geom_text(size = 3) +
    theme(axis.text.x = element_blank())
  
  
  # legend.pos <- list(
  #   xmin = 900, xmax = 1200,
  #   ymin = -1000, ymax = -1250
  # )
  # combine all plots into one
  gene.plots |>
    reduce(`+`, .init = p.path) +
    with(
      legend.pos,
      annotation_custom(
        ggplotGrob(p.legend),
        xmin = xmin, xmax = xmax,
        ymin = ymin, ymax = ymax
      )
    )
}
  
################################################################################
# https://www.kegg.jp/kegg-bin/show_pathway?syp01230
# Global amino acid path
  
foo <- my.extra.path(
  'syp01230',
  legend.pos = list(
    xmin = 900, xmax = 1200,
    ymin = -1000, ymax = -1250
  )
) +
  # remove ugly thick border
  theme(plot.margin = unit(c(-1, -2, -1, -2), "cm"))

ggsave('foo.jpeg', plot = foo,
       width = 12, height = 12, dpi = 400)

########################################
################################################################################
################################################################################