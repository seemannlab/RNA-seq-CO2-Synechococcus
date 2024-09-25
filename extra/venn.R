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
    '0.04->4% CO2', '4->30% CO2'
  )] |>
  set_names(
    '0.04% vs 4%',
    '4% vs 30%'
  ) -> de.list

'analysis/K_logFC-vs-30.tsv' |>
  read_tsv() |>
  filter(abs(log2FoldChange) >= 1, padj <= 0.001) |>
  # filter(padj <= 0.05) |>
  pull(Geneid) -> de.list$`30% vs rest`

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
  ggtitle('|log2FC| ≥ 1')

ggsave('extra/venn.jpeg', width = 8, height = 8)


################################################################################

xs <- setdiff(
  de.list$`30% vs rest`,
  c(de.list$`0.04% vs 4%`, de.list$`4% vs 30%`)
)

'data/C_annotation.tsv' |>
  read_tsv() |>
  filter(Geneid %in% xs) |>
  transmute(
    Geneid,
    nice = sprintf(
      '%s (%s)',
      product, 
      old_locus_tag |> str_remove('SYNPCC7002_')
    )
  ) |>
  arrange(nice) -> foo

write_tsv(foo, 'extra/venn-extra.tsv')

'analysis/D_normalized-counts.tsv' |>
  read_tsv() |>
  right_join(foo, 'Geneid') |>
  select(- Geneid) |>
  pivot_longer(- nice, names_to = 'lib') |>
  left_join(
    'data/C_meta.tsv' |> read_tsv(),
    'lib'
  ) -> bar

bar |>
  mutate_at('nice', str_replace, '(?<=.{20}) ', '\n') |>
  mutate_at('nice', str_replace, '(?<=[^\n]{20}) ', '\n') |>
  mutate(CO2 = CO2 |> as.character() |>fct_reorder(CO2)) |>
  ggplot(aes(as.factor(CO2), value, color = CO2)) +
  geom_jitter(width = .1) +
  xlab(NULL) +
  # scale_y_log10() +
  ylab('Expression count, DESeq2 normalized') +
  ggsci::scale_color_jama() +
  facet_wrap(~ nice, scales = 'free_y') +
  theme_pubr(12) +
  theme(legend.position = 'hide')

ggsave('extra/venn-extra.jpeg', width = 16, height = 9)
