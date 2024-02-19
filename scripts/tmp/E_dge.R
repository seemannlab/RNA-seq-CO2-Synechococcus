#!/usr/bin/env Rscript

# Differential expression analysis
library(tidyverse)
library(ggpubr)
library(patchwork)

library(DESeq2)
library(corrplot)

library(conflicted)
conflicts_prefer(dplyr::filter)
conflicts_prefer(dplyr::rename)
conflicts_prefer(purrr::reduce)
conflicts_prefer(dplyr::count)

ALPHA <- 0.05

# https://riptutorial.com/r/example/28354/colorblind-friendly-palettes
cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442",
               "#0072B2", "#D55E00", "#CC79A7")

################################################################################
# Load input data

meta <- read_tsv('data/C_meta.tsv') %>%
  mutate_at('condition', ~ fct_reorder(as.character(.x), as.numeric(.x)))
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
# DESeq2 stage-wise testing for differential expression of all pair-wise
# comparisons

DESeq2::DESeqDataSetFromMatrix(
  countData = raw.noribo.mat,
  colData = meta,
  # ! it is not possible to add a the 'batch' here, because otherwise
  #   the design matrix would not have a full rank.
  design = ~ condition
) -> des

# Screen for changing expression -> LRT
des %>%
  DESeq2::DESeq(
    test = 'LRT',
    reduced = ~ 1
  ) %>%
  DESeq2::results(alpha = ALPHA, tidy = TRUE) -> des.screen

# Confirm which pair changes -> WALD
des.confirm <- DESeq2::DESeq(des)

# Save normalized expression
des.norm <- DESeq2::counts(des.confirm, normalized = TRUE)
des.norm %>%
  as_tibble(rownames = 'Geneid') -> des.norm.tbl
write_tsv(des.norm.tbl, 'analysis/E_normalized-counts.tsv')

# Compute P-values for each pair-wise comparison
# List all pairwise comparison in format 'x->y' with
# 1. x < y
# 2. logFC > 0 means higher expression in y
meta %>%
  pull(condition) %>%
  levels %>%
  crossing(x = ., y = .) %>%
  filter(as.numeric(x) < as.numeric(y)) %>%
  with(set_names(
    # the condition it changes to `y` is the numerator, so before `x` in the contrast
    map2(x, y, ~ c('condition', .y, .x)),
    sprintf('%s->%s%% CO2', x, y)
  )) %>%
  # now match to P-Values and logFCs
  map(~ DESeq2::results(
    contrast = .x,
    object = des.confirm,
    alpha = ALPHA,
    tidy = TRUE
  )) %>%
  map2(names(.), ~ mutate(.x, test = .y)) %>%
  bind_rows() %>%
  as_tibble() -> confirm.tbl

# # Visually confirm that the logFC are in the 'right' direction
# # positive logFC -> larger expression in higher CO2 condition
# confirm.tbl %>%
#   drop_na(log2FoldChange) %>%
#   arrange(log2FoldChange) %>%
#   tail
# 
# des.norm.tbl %>%
#   # filter(Geneid == 'gene-SYNPCC7002_RS12370') %>%
#   filter(Geneid == 'gene-SYNPCC7002_RS12370') %>%
#   gather('lib', 'v', - Geneid) %>%
#   left_join(meta, 'lib') %>%
#   ggplot(aes(condition, v)) +
#   geom_boxplot() +
#   scale_y_log10()

# Ensure that P-values are well behaved
# confirm.tbl %>%
#   ggplot(aes(pvalue)) +
#   geom_histogram() +
#   facet_wrap(~ test)

################################################################################
# Stage-wise correction
# Use FDR on screen and family-wise on confirm to guarantee overall FDR

# Screen: LRT as vector
des.screen %>%
  with(set_names(
    padj, row
  )) -> screen.vec

# Confirm: Pvalues of pairwise as matrix
confirm.tbl %>%
  select(row, test, pvalue) %>%
  spread(test, pvalue) -> foo
foo %>%
  select(- row) %>%
  as.matrix %>%
  magrittr::set_rownames(foo$row) -> confirm.mat

# screen/confirm in same order
confirm.mat <- confirm.mat[names(screen.vec), ]

# Stage-wise testing
multi.test <- stageR::stageR(pScreen = screen.vec,
                             pConfirmation = confirm.mat,
                             pScreenAdjusted = TRUE)

multi.adj <- stageR::stageWiseAdjustment(multi.test,
                                         method = 'holm',
                                         allowNA = TRUE,
                                         alpha = ALPHA)

# Prepare tidy table
stageR::getAdjustedPValues(multi.adj,
                           order = FALSE,
                           onlySignificantGenes = FALSE) -> ps.stage
ps.stage %>%
  as_tibble() %>%
  mutate(row = row.names(ps.stage)) %>%
  select(- padjScreen) %>%
  gather('test', 'padj.stageR', - row) %>%
  drop_na -> ps.stage

confirm.tbl %>%
  select(row, baseMean, test, log2FoldChange, pvalue, padj) %>%
  left_join(ps.stage, c('row', 'test')) %>%
  mutate(is.de = padj.stageR <= ALPHA) %>%
  rename(Geneid = row) %>%
  left_join(annot, 'Geneid') -> deg.res

write_tsv(deg.res, 'analysis/E_dge-stagewise-analysis.tsv')

################################################################################
# Custom PCA plotting helper (same as in D_verify.R)

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
# Show relation ships in PCA

ps <- list(
  my.pca(des.confirm, 'PC1', 'PC2'),
  my.pca(des.confirm, 'PC1', 'PC3'),
  my.pca(des.confirm, 'PC1', 'PC4'),
  my.pca(des.confirm, 'PC2', 'PC3'),
  my.pca(des.confirm, 'PC2', 'PC4'),
  my.pca(des.confirm, 'PC3', 'PC4')
)

p.ps <-
  ps |>
  map(`+`, theme(legend.position = 'hide')) |>
  reduce(`+`)

p.leg <-
  cowplot::plot_grid(
    get_legend(ps[[1]]),
    ncol = 1, rel_heights = c(.9, .1)
  )

(p.ps / p.leg) +
  plot_layout(heights = c(6, 1))


ggsave('analysis/E_pca.jpeg',
       width = 18, height = 10,
       dpi = 400)


################################################################################
# Volcano plots

deg.res %>%
  ggplot(aes(log2FoldChange, - log10(padj.stageR))) +
  geom_point(alpha = 0.7) +
  geom_hline(yintercept = - log10(ALPHA), color = 'red') +
  ylab('-log10(P-value, stage-wise adjusted)') +
  facet_wrap(~ test) +
  theme_bw(18) +
  ggtitle('Volcano plots',
          paste(
            'Positive logFC is higher expression in higher CO2 concentration',
            sprintf('Red line: FDR %s', ALPHA),
            sep = '\n'
          ))

ggsave('analysis/E_volcano.jpeg',
       width = 18, height = 14,
       dpi = 400)

################################################################################
# Check for consistencies of logFCs

deg.res %>%
  select(Geneid, test, log2FoldChange) %>%
  mutate_at('log2FoldChange', replace_na, 0) %>%
  mutate_at('test', str_remove, '% CO2') %>%
  separate(test, c('from', 'to'), sep = '->') -> foo

foo %>%
  left_join(foo, c('Geneid', 'to' = 'from')) %>%
  transmute(
    Geneid, from, to = to.y,
    lfc.indirect = log2FoldChange.x + log2FoldChange.y
  ) %>%
  left_join(foo, c('Geneid', 'from', 'to')) %>%
  ggscatter(
    x = 'log2FoldChange',
    y = 'lfc.indirect',
    add = 'reg.line',
    add.params = list(color = 'blue'),
    cor.coef = TRUE,
    cor.coeff.args = list(method = "pearson", label.x = -5, label.sep = "\n")
  ) +
  xlab('"Direct" logFC') +
  ylab('"Indirect" logFC') +
  ggtitle(
    'Comparison logFC relative to third condition',
    'A->B vs A->C + C->B'
  ) +
 theme_pubr(16)

ggsave('analysis/E_indirect-logFC.jpeg',
     width = 8, height = 8,
     dpi = 400)

###############################################################################
# Overview table

foo <- 
  deg.res |>
  filter(is.de)

bar <-
  foo |>
  mutate(x = 'DE in any comparison') %>%
  bind_rows(
    . |>
      filter(! str_detect(test, '^0\\.04')) |>
      mutate(x = 'Excluding air')
  ) |>
  select(x, Geneid, type) |>
  unique() |>
  count(x, type)

foo %>%
  count(type, test) %>%
  spread(type, n) |>
  mutate_if(is.numeric, replace_na, 0) %>%
  bind_rows(
    bar |>
      rename(test = x) |>
      spread(type, n)
  ) |>
  write_tsv('analysis/E_overview-by-type.tsv')

################################################################################
# No DEGs above logFC

list(0, 1, 2, 3, 4, 5) |>
  map(function(i) {
    deg.res |>
      filter(is.de, abs(log2FoldChange) >= i) |>
      mutate(x = sprintf('|logFC| ≥ %g', i))
  }) |>
  bind_rows() -> foo
  
foo |>
  count(x, test) |>
  spread(x, n, fill = 0) |>
  bind_rows(
    foo |>
      mutate(test = 'DE in any comparison') |>
      select(test, x, Geneid) |>
      unique() |>
      count(test, x) |>
      spread(x, n, fill = 0),
    foo |>
      filter(! str_detect(test, '^0\\.04')) |>
      mutate(test = 'Excluding air') |>
      select(test, x, Geneid) |>
      unique() |>
      count(test, x) |>
      spread(x, n, fill = 0),
  ) |>
  write_tsv('analysis/E_overview-by-logFC.tsv')

################################################################################
# Heatmap of expression

deg.res %>%
  filter(is.de) %>%
  pull(Geneid) %>%
  unique -> mask

cl <- list(
  'CO2' = meta %>%
    pull(condition) %>%
    levels %>%
    set_names(cbPalette[2:(length(.) + 1)], .)
)

with(
  meta,
  data.frame('CO2' = as.character(condition), row.names = lib)
) -> cl.df

norm.vst <-
  des.confirm |>
  DESeq2::vst() |>
  SummarizedExperiment::assay()

norm.vst |>
  as_tibble(rownames = 'Geneid') |>
  write_tsv('analysis/E_vst.tsv')

# Cluster columns ahead of pheatmap to rotate tree nicely
dat <-  norm.vst[mask, ]

lib.clust <-
  dat |>
  t() |>
  dist() |>
  hclust() |>
  ape::as.phylo() |>
  ape::rotateConstr(
    meta |>
      arrange(condition, sample) |>
      pull(lib)
  ) |>
  as.hclust()

pheatmap::pheatmap(
  dat,
  scale = 'row',
  cluster_cols = lib.clust,
  show_rownames = FALSE,
  show_colnames = FALSE,
  annotation_col = cl.df,
  annotation_colors = cl,
  color = colorRampPalette(rev(
    RColorBrewer::brewer.pal(n = 7, name = "RdBu")))(59),
  filename = 'analysis/E_heatmap.jpeg'
)

################################################################################
# Correlation of logFCs

deg.res |>
  filter(is.de) |>
  select(Geneid) |>
  unique() |>
  left_join(deg.res, 'Geneid') |>
  select(Geneid, test, log2FoldChange) %>%
  spread(test, log2FoldChange, fill = 0) -> lfc

lfc %>%
  select(- Geneid) %>%
  as.matrix %>%
  magrittr::set_rownames(lfc$Geneid) -> lfc.mat

# logFC mat in which at least one comparison has abs logFC > 3
deg.res |>
  filter(is.de, abs(log2FoldChange) >= 3) |>
  pull(Geneid) |>
  unique() -> mask
lfc2 <- lfc.mat[mask, ]

jpeg('analysis/E_logFC-cor.jpeg', res = 400,
     width = 20, height = 20, units = 'cm')
lfc.mat |>
  cor() |>
  corrplot(
    method = 'square',
    order = 'original',
    type = 'lower',
    diag = FALSE,
    insig='blank',
    addCoef.col ='black',
    number.cex = 0.8,
    col = colorRampPalette(c('#FF6666', '#DDDDDD', '#6666FF'))(50)
  )
dev.off()

jpeg('analysis/E_logFC-min3-cor.jpeg', res = 400,
     width = 20, height = 20, units = 'cm')
lfc2 |>
  cor() |>
  corrplot(
    method = 'square',
    order = 'original',
    type = 'lower',
    diag = FALSE,
    insig='blank',
    addCoef.col ='black',
    number.cex = 0.8,
    col = colorRampPalette(c('#FF6666', '#DDDDDD', '#6666FF'))(50)
  )
dev.off()

################################################################################
# Heatmap all logFCs


helper <- function(lfc, path) {
  myColors <- colorRampPalette(rev(RColorBrewer::brewer.pal(n = 7, name = "RdBu")))(59)
  myBreaks <- c(seq(min(lfc), 0, length.out = ceiling(length(myColors) / 2) + 1), 
                seq(max(lfc) / length(myColors), max(lfc),
                    length.out = floor(length(myColors) / 2)))
  
  pheatmap::pheatmap(
    lfc,
    scale = 'none',
    show_rownames = FALSE,
    color = myColors,
    breaks = myBreaks,
    filename = path
  )
}

helper(lfc.mat, 'analysis/E_logFC-heatmap.jpeg')
helper(lfc2, 'analysis/E_logFC-min3-heatmap.jpeg')
