#!/usr/bin/env Rscript

# Size-factor norm the public data and show the general structure
library(tidyverse)
library(ggpubr)
library(patchwork)

################################################################################
# Load input data

public.meta <-
  'data/C_public-meta.tsv' |>
  read_tsv() |>
  left_join(
    'data/A_public-stat.tsv' |>
      read_tsv() |>
      select(platform = Platform, sample = SRA),
    'sample'
  )

annot <- read_tsv('data/C_annotation.tsv')
public.raw <- read_tsv('data/C_public-raw-counts.tsv')

public.raw %>%
  select(- Geneid) %>%
  as.matrix() %>%
  magrittr::set_rownames(public.raw$Geneid) -> raw.mat


################################################################################
# Exclude rRNA

annot %>%
  filter(type != 'rRNA') %>%
  pull(Geneid) -> mask

raw.noribo <- raw.mat[mask, ]

################################################################################
# DESeq2 size-factor norm and vst

des <-
  DESeq2::DESeqDataSetFromMatrix(
    countData = raw.noribo,
    colData = public.meta,
    # ! it is not possible to add a the 'batch' here, because otherwise
    #   the design matrix would not have a full rank.
    design = ~ I(paste(condition, dataset))
  ) |>
  DESeq2::DESeq()

# SF normalized expression
des.norm <- DESeq2::counts(des, normalized = TRUE)
des.norm %>%
  as_tibble(rownames = 'Geneid') |>
  write_tsv('analysis/F_public-normalized-counts.tsv')

des.vst <-
  des |>
  DESeq2::vst() |>
  SummarizedExperiment::assay()

des.vst |>
  as_tibble(rownames = 'Geneid') |>
  write_tsv('analysis/F_public-vst.tsv')

################################################################################
# PCA plot to get an idea of the overall structure

# Adapted from DESeq2::plotPCA
ntop <- 500
rv <- MatrixGenerics::rowVars(des.vst)
sel <- order(rv, decreasing = TRUE)[seq_len(min(ntop, length(rv)))]
pca <- prcomp(t(des.vst[sel, ]))
percentVar <- ( pca$sdev^2/sum(pca$sdev^2) ) %>%
  set_names(paste0('PC', 1:length(.)))

x <- 'PC1'
y <- 'PC2'
pca$x |>
  as_tibble(rownames = 'lib') |>
  left_join(public.meta, 'lib') |>
  ggplot(aes(!! sym(x), !! sym(y),
             shape = platform, color = dataset)) +
  geom_point(size = 5) +
  xlab(sprintf('%s: %.1f%% variance', x, percentVar[[x]] * 100)) +
  ylab(sprintf('%s: %.1f%% variance', y, percentVar[[y]] * 100)) +
  ggsci::scale_color_jco() +
  coord_cartesian() +
  guides(color = guide_legend(nrow = 3)) +
  theme_pubr(18) +
  theme(legend.box = 'vertical')

 ggsave('analysis/F_public-pca.jpeg', width = 13, height = 8)
