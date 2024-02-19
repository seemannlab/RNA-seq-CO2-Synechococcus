#!/usr/bin/env Rscript

# Make the high/low comparisons
# That is test the 1->15% and 4/8->30% comparison

library(tidyverse)
library(ggpubr)
library(patchwork)

library(venn)
library(corrplot)

library(conflicted)
conflicts_prefer(dplyr::filter)
conflicts_prefer(dplyr::rename)
conflicts_prefer(purrr::reduce)
conflicts_prefer(dplyr::count)

source('scripts/helper_deg.R')

ALPHA <- 0.05

# https://riptutorial.com/r/example/28354/colorblind-friendly-palettes
cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442",
               "#0072B2", "#D55E00", "#CC79A7")

################################################################################
# Load input data

meta <- read_tsv('data/C_meta.tsv') %>%
  mutate_at(c('CO2', 'photons'), ~ fct_reorder(as.character(.x), as.numeric(.x)))
annot <- read_tsv('data/C_annotation.tsv')
raw.counts <- read_tsv('data/C_raw-counts.tsv')

raw.counts %>%
  select(- Geneid) %>%
  as.matrix() %>%
  magrittr::set_rownames(raw.counts$Geneid) -> raw.counts.mat


################################################################################
# Exclude rRNA

annot %>%
  filter(type != 'rRNA') %>%
  pull(Geneid) -> mask

raw.noribo.mat <- raw.counts.mat[mask, ]

################################################################################
# Build DESeq2 Datasets per light intensity batch

dat <-
  meta |>
  pull(photons) |>
  unique() %>%
  set_names(as.character(.)) |>
  map(~ filter(meta, photons == .x)) |>
  map(function(x) {
    DESeq2::DESeqDataSetFromMatrix(
      countData = raw.noribo.mat[, x$lib],
      colData = x,
      design = ~ CO2
    )
  })

################################################################################

des <- 
  dat |>
  map(DESeq2::DESeq)

x <- des$`250` |> DESeq2::vst()
p1 <- helper.pca(x, sh = 'batch', color_values = cbPalette[-5][- c(2, 5)])

x2 <- des$`150` |> DESeq2::vst()
p2 <- helper.pca(x2, sh = 'batch', color_values = cbPalette[-5][c(2, 5)])

################################################################################

p1 + p2 +
  plot_annotation(title = 'PCA on vst data per dataset (analysed separately)')

ggsave('analysis/D_PCA-per-dataset.jpeg', width = 18, height = 8)

################################################################################

todo <- list(
  '1 -> 15' = list(photons = '150', contrast = c('CO2', 15, 1)),
  '4 -> 30' = list(photons = '250', contrast = c('CO2', 30, 4)),
  '8 -> 30' = list(photons = '250', contrast = c('CO2', 30, 8)),
  'average of 4/8 -> 30' = list(photons = '250', contrast = c(
    "Intercept"      = 0,
    "CO2_4_vs_0.04"  = - 0.5,
    "CO2_8_vs_0.04"  = - 0.5,
    "CO2_30_vs_0.04" = 1
  ))
)


high.low <-
  todo %>%
  map2(names(.), function(x, xn) {
    des[[x$photons]] |>
      DESeq2::results(
        contrast = x$contrast,
        alpha = ALPHA,
        tidy = TRUE
      ) |>
      as_tibble() |>
      mutate(cmp = xn)
  }) |>
  bind_rows()

write_tsv(high.low, 'analysis/D_dge-per-dataset.tsv')
  
################################################################################

high.low |>
  select(row, cmp, log2FoldChange) |>
  spread(cmp, log2FoldChange) |> 
  select(- row) |>
  GGally::ggpairs() +
  theme_bw(18)

ggsave('analysis/D_correlation-logFCs.jpeg', width = 12, height = 12)

################################################################################

high.low |>
  filter(padj <= ALPHA) |>
  filter(cmp %in% c('1 -> 15', 'average of 4/8 -> 30')) |>
  mutate(cmp = paste(cmp, ifelse(
    log2FoldChange > 0,
    'Up-reg.',
    'Down-reg.'
  ), sep = '\n')) |>
  group_by(cmp) |>
  do(i = list(.$row)) |>
  with(set_names(i, cmp)) |>
  map(1) |>
  venn(zcolor = 'style', ilcs = 1.5, sncs = 1.5,
       box = FALSE, ggplot = TRUE) +
  annotate('text', 50, 950,
           # label = paste('FDR', j, '%, abs logFC ≥', i),
           label = 'FDR 5%',
           hjust = 0, size = 8) -> p1

high.low |>
  filter(padj <= ALPHA, abs(log2FoldChange) >= 1) |>
  filter(cmp %in% c('1 -> 15', 'average of 4/8 -> 30')) |>
  mutate(cmp = paste(cmp, ifelse(
    log2FoldChange > 0,
    'Up-reg.',
    'Down-reg.'
  ), sep = '\n')) |>
  group_by(cmp) |>
  do(i = list(.$row)) |>
  with(set_names(i, cmp)) |>
  map(1) |>
  venn(zcolor = 'style', ilcs = 1.5, sncs = 1.5,
       box = FALSE, ggplot = TRUE) +
  annotate('text', 50, 950,
           label = 'FDR 5% & |logFC| ≥ 1',
           hjust = 0, size = 8) -> p2


p1 + p2
ggsave('analysis/D_venn-per-dataset.jpeg', width = 20, height = 10)


################################################################################
################################################################################

################################################################################
################################################################################