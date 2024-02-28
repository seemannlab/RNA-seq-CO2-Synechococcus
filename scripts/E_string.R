# Export data for use in stringApp

library(tidyverse)
library(ggpubr)

################################################################################
# Load data

annot <-
  'data/C_annotation.tsv' |>
  read_tsv()

dge <-
  'analysis/D_stagewise-adjusted-DEGs.tsv' |>
  read_tsv()

vst <-
  'analysis/D_vst-expression.tsv' |>
  read_tsv()

meta <-
  'data/C_meta.tsv' |>
  read_tsv()


################################################################################
# Loci of DE genes

de.loci <-
  dge |>
  filter(is.de) |>
  select(Geneid, locus = old_locus_tag) |>
  drop_na() |>
  unique()

################################################################################
# Z-scale VST expression with subsequent averaging

vst.de <-
  vst |>
  semi_join(
    dge |>
      filter(is.de) |>
      select(Geneid) |>
      unique(),
    'Geneid'
  )

vst.de.mat <-
  vst.de |>
  select(- Geneid) |>
  as.matrix() |>
  magrittr::set_rownames(vst.de$Geneid)

z.mat <-
  vst.de.mat |>
  apply(1, scale) |>
  t() |> 
  magrittr::set_colnames(colnames(vst.de.mat))

z.expr <-
  z.mat |>
  as_tibble(rownames = 'Geneid') |>
  pivot_longer(- Geneid, names_to = 'lib', values_to = 'z.expression') |>
  left_join(meta, 'lib')

z.avg <-
  z.expr |>
  group_by(Geneid, CO2) |>
  summarize(avg.z = mean(z.expression)) |>
  ungroup()

################################################################################
# Determine a cutoff for association scores

string <-
  'https://stringdb-downloads.org/download/protein.links.detailed.v12.0/32049.protein.links.detailed.v12.0.txt.gz' |>
  read_delim(delim = ' ') |>
  mutate_if(is.numeric, ~ .x / 1e3)

string |>
  pull(combined_score) |>
  quantile(.8) -> xs

string |>
  ggplot(aes(combined_score)) +
  scale_x_continuous(breaks = seq(0, 1, .1)) +
  scale_y_continuous(breaks = seq(0, 1, .1)) +
  geom_hline(yintercept = .8, color = 'red') +
  geom_vline(xintercept = xs, color = 'red') +
  annotate('text', x = .75, y = .1, label = xs, color = 'red', size = 5) +
  stat_ecdf() +
  xlab('STRING association score') +
  ylab('Empirical cumulative density') +
  theme_pubr(18)

ggsave('analysis/E_string-cutoff.jpeg', width = 7, height = 6, dpi = 400)

################################################################################
# Export

de.loci |>
  pull(locus) |>
  unique() |>
  write_lines('analysis/E_string-loci.txt')

z.avg |>
  inner_join(de.loci, 'Geneid') |>
  mutate(locus = paste0('32049.', locus)) |>
  write_tsv('analysis/E_string-z-expression.tsv')


################################################################################
sessionInfo()