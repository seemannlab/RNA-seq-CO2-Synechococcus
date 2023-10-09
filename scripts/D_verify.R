#!/usr/bin/env Rscript

# RNA concentrations vary significantly
# -> verify DESeq2 normalization

library(tidyverse)
library(ggpubr)
library(patchwork)

library(UpSetR)

library(furrr)
plan(multisession)

cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442",
               "#0072B2", "#D55E00", "#CC79A7")


ALPHA <- 0.05

################################################################################
# Genesets for genes up to a varying stability level

xs.pct <- c(1, 5, 10, 20, 50, 100)

xs.lab <- sprintf('Top %g%%', xs.pct)
xs.lab[length(xs.lab)] <- 'All genes'

data('dat.stability', package =  'GenomescaleStableGenes')

xs <- quantile(dat.stability$fluctuation, probs = xs.pct / 100)

xs %>%
  map(~ filter(dat.stability, fluctuation <= .x)) %>%
  map(pull, 'Geneid') -> ref.genes

tibble(
  'Stable coding genes' = xs.lab, 
  'Fluctuation cutoff' = unlist(xs),
  'No. genes' = map(ref.genes, length) %>% unlist
) %>%
  write_tsv('analysis/D_refgenes-overview.tsv')

################################################################################
# Project data

raw.counts <-
  'data/C_raw-counts.tsv' |>
  read_tsv()

meta <- 
  'data/C_meta.tsv' |>
  read_tsv() |>
  mutate_at('condition', as.factor)

################################################################################
# Expression matrix for coding genes only

raw.counts %>%
  select(- Geneid) %>%
  as.matrix %>%
  magrittr::set_rownames(raw.counts$Geneid) -> raw.counts.mat

coding.mat <- raw.counts.mat[dat.stability$Geneid, ]

################################################################################
# For different reference genes, estimate size-factors etc

custom.norm <- function(refs, libs = NULL) {
  # Potential subset by libraries (eg for excluding 0.04 condition)
  if (is.null(libs)) {
    sub.meta <- meta
    sub.coding <- coding.mat
  } else {
    sub.meta <- filter(meta, lib %in% libs)
    sub.coding <- coding.mat[, sub.meta$lib]
  }
  
  # Subset of reference genes
  sub.mat <- sub.coding[refs, ]
  
  # Estimate size-factors on subset of reference genes
  DESeq2::DESeqDataSetFromMatrix(
    countData = sub.mat,
    colData = sub.meta,
    design = ~ condition
  ) %>%
    DESeq2::DESeq() -> sub.des
  
  
  # Normalize dataset with size-factors of subset
  DESeq2::DESeqDataSetFromMatrix(
    countData = sub.coding,
    colData = sub.meta,
    design = ~ condition
  ) -> des.lrt
  DESeq2::sizeFactors(des.lrt) <- DESeq2::sizeFactors(sub.des)
  des.lrt <- DESeq2::DESeq(des.lrt, reduced = ~ 1, test = 'LRT')
  
  
  list(
    sizeFactors = DESeq2::sizeFactors(sub.des),
    des.lrt = des.lrt,
    des.lrt.res =  DESeq2::results(des.lrt, tidy = TRUE, alpha = ALPHA)
  )
}

ref.genes %>%
  future_map(custom.norm, .options = furrr_options(seed = 123)) %>%
  set_names(xs.lab) -> ref.norm

# same but without the air condition
noair.libs <-
  meta |>
  filter(condition != '0.04') |>
  pull(lib)

ref.genes %>%
  future_map(custom.norm,
             libs = noair.libs,
             .options = furrr_options(seed = 123)) %>%
  set_names(xs.lab) -> ref.norm.noair

################################################################################
# Compare normalization and impact on 
# - size factors
# - PCA
# - differential expression

ref.norm %>%
  map('sizeFactors') %>%
  map2(names(.), ~ tibble(lib = names(.x), sf = .x, norm = .y)) %>%
  bind_rows() %>%
  mutate_at('norm', as_factor) %>%
  left_join(meta, 'lib') %>%
  mutate(lab = ifelse(norm == 'All genes', sample, NA_character_)) %>%
  ggplot(aes(norm, sf, group = lib, color = condition, label = lab)) +
  geom_line(size = 1.4, alpha = .7) +
  ggrepel::geom_text_repel(direction = 'y', nudge_x = .5, show.legend = FALSE) +
  scale_color_manual(name = 'CO2 %',
                     values = cbPalette[-5]) +
  xlab('Normalized to gene stability subset') +
  ylab('DESeq2 size-factor') +
  theme_pubr(16) -> p1

ref.norm.noair %>%
  map('sizeFactors') %>%
  map2(names(.), ~ tibble(lib = names(.x), sf = .x, norm = .y)) %>%
  bind_rows() %>%
  mutate_at('norm', as_factor) %>%
  left_join(meta, 'lib') %>%
  mutate(lab = ifelse(norm == 'All genes', sample, NA_character_)) %>%
  ggplot(aes(norm, sf, group = lib, color = condition, label = lab)) +
  geom_line(size = 1.4, alpha = .7) +
  ggrepel::geom_text_repel(direction = 'y', nudge_x = .5, show.legend = FALSE) +
  scale_color_manual(name = 'CO2 %',
                     values = cbPalette[-c(1, 5)]) +
  xlab('Normalized to gene stability subset') +
  ylab('DESeq2 size-factor') +
  theme_pubr(16) -> p2

p1 + p2 +
  plot_annotation(tag_levels = 'A')

ggsave('analysis/D_size-factors.jpeg', width = 20, height = 10, dpi = 400)

################################################################################
# Custom PCA plotting helper

my.pca <- function(dat, x = 'PC1', y = 'PC2',
                   color_values = cbPalette[-5], ntop = 500) {
  # Adapted from DESeq2::plotPCA
  des.blog <- DESeq2::vst(dat, blind = TRUE)
  dba <- SummarizedExperiment::assay(des.blog)
  rv <- MatrixGenerics::rowVars(dba)
  sel <- order(rv, decreasing = TRUE)[seq_len(min(ntop, length(rv)))]
  pca <- prcomp(t(dba[sel, ]))
  percentVar <- ( pca$sdev^2/sum(pca$sdev^2) ) %>%
    set_names(paste0('PC', 1:length(.)))
  
  pca$x |>
    as_tibble(rownames = 'lib') |>
    left_join(meta, 'lib') |>
    ggplot(aes(!! sym(x), !! sym(y), shape = batch, color = condition)) +
    geom_point(size = 5) +
    xlab(sprintf('%s: %.1f%% variance', x, percentVar[[x]] * 100)) +
    ylab(sprintf('%s: %.1f%% variance', y, percentVar[[y]] * 100)) +
    ggrepel::geom_label_repel(aes(label = sample), show.legend = FALSE) +
    scale_color_manual(name = 'CO2 %',
                       values = color_values) +
    coord_cartesian() +
    theme_pubr(18)
}

################################################################################

# PCA for full dataset

ref.norm %>%
  map('des.lrt') %>%
  map(my.pca, x = 'PC1', y = 'PC2',
      color_values = cbPalette[-c(5)]) %>%
  map2(names(.), ~ .x + ggtitle(.y)) -> ps

ps %>%
  map(~ .x + theme(legend.position = 'hide')) %>%
  invoke(.f = cowplot::plot_grid) %>%
  cowplot::plot_grid(
    get_legend(ps[[1]]),
    ncol = 1, rel_heights = c(.9, .1)
  )

ggsave('analysis/D_pca-norm-coding-genes-with-air.jpeg',
       width = 18, height = 12, dpi = 400)

# PCA for dataset without out air

ref.norm.noair %>%
  map('des.lrt') %>%
  map(my.pca, x = 'PC1', y = 'PC2',
      color_values = cbPalette[-c(1, 5)]) %>%
  map2(names(.), ~ .x + ggtitle(.y)) -> ps

ps %>%
  map(~ .x + theme(legend.position = 'hide')) %>%
  invoke(.f = cowplot::plot_grid) %>%
  cowplot::plot_grid(
    get_legend(ps[[1]]),
    ncol = 1, rel_heights = c(.9, .1)
  )

ggsave('analysis/D_pca-norm-coding-genes-without-air.jpeg',
       width = 18, height = 12, dpi = 400)

################################################################################

ref.norm %>%
  map('des.lrt.res') %>%
  map2(names(.), ~ mutate(.x, norm = .y)) %>%
  bind_rows() -> res

res %>%
  transmute(row, norm, de = as.integer(padj <= ALPHA) %>% replace_na(0)) %>%
  spread(norm, de) -> foo
foo %>%
  select(- row) %>%
  as.data.frame() %>%
  magrittr::set_rownames(foo$row) -> bar


jpeg('analysis/D_upset.jpeg', width = 10, height = 8, res = 400, units = 'in')
upset(bar, nsets = ncol(bar), 
      nintersects = NA, order.by = 'freq',
      point.size = 3, line.size = 1.2,
      text.scale =  1.2)
dev.off()


################################################################################