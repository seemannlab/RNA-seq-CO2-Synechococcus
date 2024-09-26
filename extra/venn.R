# Illustrate the issue of up/down regulation discrepancies in comparison
# to what the 30% vs rest analysis shows
library(tidyverse)
library(patchwork)
library(ggpubr)

################################################################################
'analysis/D_stagewise-adjusted-DEGs.tsv' |>
  read_tsv() |>
  filter(is.de) |>
  filter(abs(log2FoldChange) >= 1) |>
  group_by(test) |>
  group_split() |>
  as.list() ->foo
set_names(
  foo |> map(pull, Geneid),
  foo |> map(pull, test) |> map(first) |> unlist()
) |>
  _[c(
    '0.04->4% CO2', '4->30% CO2',
    '0.04->30% CO2'
  )] |>
  set_names(
    '0.04% vs 4%',
    '4% vs 30%',
    '0.04% vs 30%',
  ) -> de.list

de.list |>
  venn::venn(
    ilabels = 'counts',
    ilcs = 2, sncs = 2,
    zcolor = 'style',
    box = FALSE,
    ggplot = TRUE
  ) -> p

wrap_elements(full = p) +
  theme_void(18) +
  ggtitle('|log2FC| ≥ 1, FDR 0.05')

ggsave('extra/venn-pairs.jpeg', width = 8, height = 8)

################################################################################

'analysis/K_logFC-vs-30.tsv' |>
  read_tsv() |>
  filter(abs(log2FoldChange) >= 1, padj <= 0.05) |>
  pull(Geneid) -> de.list$`30% vs\n0.04%,\n4%, 8%`

de.list |>
  venn::venn(
    ilabels = 'counts',
    ilcs = 2, sncs = 2,
    zcolor = 'style',
    box = FALSE,
    ggplot = TRUE
  ) -> p

################################################################################

wrap_elements(full = p) +
  theme_void(18) +
  ggtitle('|log2FC| ≥ 1, FDR 0.05')

ggsave('extra/venn-30extra.jpeg', width = 8, height = 8)


