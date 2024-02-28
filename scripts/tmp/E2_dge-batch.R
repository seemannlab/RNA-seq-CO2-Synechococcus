#!/usr/bin/env Rscript

# Repeat Differential expression analysis
# - without air
# - without first batch
# compare to full DESeq2+stageR run

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
# Helper function: Adapted from E_dge.R
# DESeq2 stage-wise testing for differential expression of all pair-wise
# comparisons

helper <- function(sub.meta, sub.mat) {
  DESeq2::DESeqDataSetFromMatrix(
    countData = sub.mat,
    colData = sub.meta,
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
  
  # Compute P-values for each pair-wise comparison
  # List all pairwise comparison in format 'x->y' with
  # 1. x < y
  # 2. logFC > 0 means higher expression in y
  sub.meta %>%
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
}

################################################################################

# subset by batch, without 1/15%
sub.meta <-
  meta |>
  filter(str_detect(batch, 'batch3')) |>
  mutate_at('condition', fct_drop)

deg.batch <- helper(sub.meta, raw.noribo.mat[, sub.meta$lib])

# subset without air

sub.meta <-
  meta |>
  filter(condition != '0.04') |>
  mutate_at('condition', fct_drop)

deg.air <- helper(sub.meta, raw.noribo.mat[, sub.meta$lib])

################################################################################

deg.batch |>
  write_tsv('analysis/E2_stagewise-batch3-only.tsv')
deg.air |>
  write_tsv('analysis/E2_stagewise-no-air.tsv')

################################################################################

deg.full <-
  'analysis/E_dge-stagewise-analysis.tsv' |>
  read_tsv()

################################################################################
# compare logFCs overall

cmp <- list(
  'Full dataset' = deg.full,
  'Excluding 1/15%' = deg.batch,
  'Excluding air' = deg.air
)


###############################################################################
# pairwise logFC scatter plots

cmp.pairs <-
  cmp |>
  map(select, Geneid, test, log2FoldChange) %>%
  map2(names(.), ~ mutate(.x, x = .y)) |>
  bind_rows() %>%
  left_join(., ., c('Geneid', 'test'), relationship = "many-to-many") |>
  filter(x.x < x.y)

cmp.pairs |>
  mutate(
    x.x = paste('x =', x.x),
    x.y = paste('y =', x.y)
  ) |>
  ggscatter(
    'log2FoldChange.x', 'log2FoldChange.y',
    add = 'reg.line', add.params = list(color = 'blue'),
    cor.coef = TRUE, cor.coeff.args = list(color = 'red')
  ) +
  facet_wrap(x.x ~ x.y)

ggsave('analysis/E2_logFC-comparisons.jpeg',
       width = 12, height = 8, dpi = 400, scale = 0.7)

###############################################################################
# sanitize logFC of degs

cmp |>
  map(filter, is.de) |>
  map(select, Geneid, test, log2FoldChange) %>%
  map2(names(.), ~ mutate(.x, x = .y)) |>
  bind_rows() %>%
  transmute(x, Geneid, test, slfc = sign(log2FoldChange)) %>%
  left_join(., ., c('Geneid', 'test'), relationship = "many-to-many") |>
  filter(x.x < x.y) |>
  count(x.x, x.y, slfc.x == slfc.y)
# x.x             x.y           `slfc.x == slfc.y`     n
# Excluding 1/15% Excluding air TRUE                3752
# Excluding 1/15% Full dataset  TRUE               10573
# Excluding air   Full dataset  TRUE               16754
  
# Perfect, all degs are also guaranteed to have the same regulation direction

###############################################################################

cmp.diff <-
  cmp |>
  map(filter, is.de) |>
  map(select, Geneid, test) %>%
  map2(names(.), ~ mutate(.x, x = .y)) |>
  bind_rows()

cmp.diff |>
  select(- test) |>
  unique() |>
  group_by(x) |>
  do(i = list(.$Geneid)) |>
  with(set_names(i, x)) |>
  map(1) |>
  venn::venn(
    zcolor = "style",
    ggplot = TRUE,
    ilcs = 1.5,
    sncs = 1.5,
    box = FALSE
  ) +
  plot_annotation(
    'Differentially expressed in any pairwise comparison',
    theme = theme_pubr(16)
  )

ggsave('analysis/E2_venn-overall.jpg',
      width = 8, height = 8, dpi = 400)


###############################################################################
conflicts_prefer(base::intersect)

list(
  # overall number of genes
  cmp.diff |>
    mutate(x = str_remove(x, 'Excluding ')) |>
    mutate(x = str_remove(x, ' dataset')) |>
    mutate(x = str_replace(x, 'Full', 'full')) %>%
    count(test, x) |>
    spread(x, n, fill = 0),
  
  # shared in pairwise comparison
  cmp.diff |>
    mutate(x = str_remove(x, 'Excluding ')) |>
    mutate(x = str_remove(x, ' dataset')) |>
    mutate(x = str_replace(x, 'Full', 'full')) %>%
    left_join(., ., c('Geneid', 'test'), relationship = "many-to-many") |>
    filter(x.x > x.y) |>
    count(test, x.x, x.y) |>
    transmute(
      test,
      n,
      y = paste(x.x, '&', x.y)
    ) |>
    spread(y, n, fill = 0),
  
  # genes that appear in all 3 datasets
  cmp.diff |>
    count(Geneid, test) |>
    filter(n == 3) |>
    count(test, name = 'All')
) |>
  reduce(.f = left_join, by = 'test') |>
  mutate_if(is.numeric, replace_na, 0) -> venn.like

venn.like |>
  pivot_longer(- test) |>
  mutate_if(is.character, fct_inorder) |>
  mutate(
    v2 = prettyNum(value, big.mark = ','),
    cl = ifelse(value > 1000, 'black', 'white')
  ) |>
  ggplot(aes(name, test, fill = value, label = v2, color = I(cl))) +
  geom_tile() +
  scale_fill_viridis_c(name = 'No. genes') +
  guides(fill = guide_colorbar(barwidth = unit(10, 'cm'))) +
  geom_text() +
  ylab('Significant logFC') +
  xlab('Sets and intersections') +
  theme_pubr(18)

ggsave('analysis/E2_set-intersections.jpeg',
       width = 14, height = 12, dpi = 400)

################################################################################
# A nicer comparison of logFC in the full analysis

deg <- deg.full

# genes sharing combinations of logFC across pairwise comparisons
shared.genes <-
  deg |>
  filter(is.de) |>
  transmute(
    Geneid,
    test = str_remove(test, '% CO2$'),
    reg = ifelse(log2FoldChange > 0, 'logFC > 0', 'logFC < 0')
  ) |>
  spread(test, reg, fill = 'logFC not sig.') |>
  select(- Geneid) |>
  group_by_all() |>
  tally() |>
  ungroup() |>
  arrange(desc(n)) |>
  filter(n >= 10)  |>
  mutate(row = paste0('row', 1:n()) |>
           fct_inorder() |>
           fct_rev())


# custom order of comparisons per hierarchical clustering
test.order <-
  shared.genes |>
  select(contains('->')) |>
  mutate_all(fct_relevel,
            'logFC > 0', 'logFC < 0', 'logFC not sig.') |>
  mutate_all(as.integer) |>
  as.matrix() |>
  t() |>
  dist() |>
  hclust() |>
  with(labels[order])

# plot heatmap of regulation
p1 <-
  shared.genes |>
  select(- n) |>
  pivot_longer(- row) |>
  mutate_at('name', fct_relevel, test.order ) |>
  mutate_at('value', fct_relevel,
            'logFC > 0', 'logFC < 0', 'logFC not sig.') |>
  ggplot(aes(name, row, fill = value)) +
  geom_tile(color = 'black') +
  # scale_fill_manual(values = cbPalette[c(7, 3, 1)],
  scale_fill_manual(values = c('red', 'blue', 'grey'),
                    name = NULL) +
  xlab('Pairwise expression comparison') +
  ylab(NULL) +
  theme_pubr(18) +
  theme(
    axis.text.x = element_text(angle = 60, hjust = 1),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank()
  )

p2 <-
  shared.genes |>
  mutate(nice = prettyNum(n, big.mark = ',')) |>
  ggplot(aes(n, row, label = nice)) +
  geom_bar(stat = 'identity') +
  geom_text(nudge_x = 5, size = 5) +
  ylab(NULL) +
  xlab('No. genes') +
  theme_pubr(18) +
  theme(
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank()
  )
  
p1 | p2

ggsave('analysis/E2_logFC-cmp.jpeg',
       width = 10, height = 12, dpi = 400)
