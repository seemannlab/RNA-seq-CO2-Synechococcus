library(tidyverse)
library(venn)

################################################################################

res1 <-
  'analysis/D_dge-per-dataset.tsv' |>
  read_tsv()

res2 <-
  'analysis/E_stagewise-adjusted-DEGs.tsv' |>
  read_tsv()

################################################################################

venn.dat <-
  res2 |>
  filter(is.de) |>
  select(version, Geneid) |>
  mutate(version = ifelse(
    version == 'Without correction',
    'Full dataset\n(stage-wise FDR 5%)',
    'Full dataset with correction\n(stage-wise FDR 5%)'
  )) %>%
  group_by(version) |>
  do(i = .$Geneid) |>
  with(set_names(i, version))

venn.dat$`1&15% only\n(single logFC FDR 5%)` <-
  res1 |>
  filter(cmp == '1 -> 15', padj <= 0.05) |>
  pull(row)
  

################################################################################

meta <-
  'data/C_meta.tsv' |>
  read_tsv() |>
  mutate_at('batch', str_remove, '-Lane.$') |>
  filter(CO2 != 1, CO2 != 15) |>
  mutate_at(c('CO2', 'photons'), function(x) {
    x |>
      as.character() |>
      fct_reorder(as.numeric(x))
  })


annot <-
  'data/C_annotation.tsv' |>
  read_tsv()

annot %>%
  filter(type != 'rRNA') %>%
  pull(Geneid) -> mask


raw <-
  'data/C_raw-counts.tsv' |>
  read_tsv()

raw.mat <-
  raw |>
  select(- Geneid) |>
  as.matrix() |>
  magrittr::set_rownames(raw$Geneid) |>
  _[mask, meta$lib]

################################################################################
# various versions of Combat correction

des <- DESeq2::DESeqDataSetFromMatrix(
  countData = raw.mat,
  colData = meta,
  design = ~ CO2
)
source('scripts/helper_deg.R')
res3 <- helper.stagewise(des, helper.pairwise.deg(des), 0.05)

write_tsv(res3, 'analysis/E2_stagewise-without-1-15.tsv')


################################################################################

venn.dat$`Without 1&15%\n(stage-wise FDR 5%)` <-
  res3 |>
  filter(is.de) |>
  pull(Geneid) |>
  unique()


################################################################################

venn(venn.dat, zcolor = 'style',
     box = FALSE, ggplot = TRUE,
     ilcs = 1.5, sncs = 1.5)

ggsave('analysis/E2_venn-adding-1-15.jpeg', width = 10, height = 10)


################################################################################

res2 |>
  select(Geneid, type, name, product, test, version, log2FoldChange) |>
  spread(version, log2FoldChange) |>
  mutate(foo = abs(`With correction` - `Without correction`))  |>
  filter(abs(`Without correction`) > 1)  |>
  arrange(foo) |>
  filter(! str_detect(product, 'hypothetical')) |>
  View()
