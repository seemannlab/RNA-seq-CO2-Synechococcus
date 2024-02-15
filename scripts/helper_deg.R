# Helper scripts for modular analysis scripts:
# - helper.pairwise.deg
# - helper.stagewise
# - helper.pca

################################################################################
################################################################################

# DESeq differential expression for all pairwise comparisons in colData$CO2
#
# The helper function input is a DESeq DataSet, created similiar to:
# DESeq2::DESeqDataSetFromMatrix(
#   countData = raw.noribo.mat,
#   colData = meta,
#   design = ~ condition
# )
# Output: Table of logFCs and p-Values with column `test` identifying the pair
helper.pairwise.deg <- function(des, ALPHA = 0.05) {
  # Compute P-values for each pair-wise comparison
  # List all pairwise comparison in format 'x->y' with
  # 1. x < y
  # 2. logFC > 0 means higher expression in y
  des <- DESeq2::DESeq(des)
  des |>
    SummarizedExperiment::colData() |>
    as_tibble() |>
    pull(CO2) |>
    levels() %>%
    crossing(x = ., y = .) |>
    filter(as.numeric(x) < as.numeric(y)) |>
    with(set_names(
      # the condition it changes to `y` is the numerator, so before `x` in the contrast
      map2(x, y, ~ c('CO2', .y, .x)),
      sprintf('%s->%s%% CO2', x, y)
    )) |>
    # now match to P-Values and logFCs
    map(~ DESeq2::results(
      contrast = .x,
      object = des,
      alpha = ALPHA,
      tidy = TRUE
    )) %>%
    map2(names(.), ~ mutate(.x, test = .y)) |>
    bind_rows() %>%
    as_tibble()
}

################################################################################
################################################################################

# Input: DEseq Dataset and output from helper.pairwise.deg
# Output: logFC table with additionally stage-wise adjusted P-Values
helper.stagewise <- function(des, pairs, ALPHA) {
  # Stage-wise correction
  # Use FDR on screen and family-wise on confirm to guarantee overall FDR
  
  # Screen for changing expression -> LRT
  des.screen <-
    des |>
    DESeq2::DESeq(
      test = 'LRT',
      reduced = ~ 1
    ) |>
    DESeq2::results(alpha = ALPHA, tidy = TRUE)
  
  # Screen: LRT P-Value as vector
  screen.vec <-
    des.screen |>
    with(set_names(
      padj, row
    ))
  
  # Confirm: Pvalues of pairwise tests as matrix
  pairs |>
    select(row, test, pvalue) |>
    spread(test, pvalue) -> foo
  foo |>
    select(- row) |>
    as.matrix() |>
    magrittr::set_rownames(foo$row) -> pairs.mat
  
  # screen/confirm in same order
  pairs.mat <-pairs.mat[names(screen.vec), ]
  
  # Stage-wise testing
  multi.test <- stageR::stageR(
    pScreen = screen.vec,
    pConfirmation = pairs.mat,
    pScreenAdjusted = TRUE
  )
  multi.adj <- stageR::stageWiseAdjustment(
    multi.test,
    method = 'holm',
    allowNA = TRUE,
    alpha = ALPHA
  )
  
  # Prepare tidy table
  ps.stage <-
    stageR::getAdjustedPValues(
      multi.adj,
      order = FALSE,
      onlySignificantGenes = FALSE
    ) |>
    as_tibble(rownames = 'row') |>
    select(- padjScreen) |>
    gather('test', 'padj.stageR', - row) |>
    drop_na()
  
  # Combine and add the annotation for tiyier results
  pairs |>
    select(row, baseMean, test, log2FoldChange, pvalue, padj) |>
    left_join(ps.stage, c('row', 'test')) |>
    mutate(is.de = padj.stageR <= ALPHA) |>
    rename(Geneid = row) |>
    left_join(annot, 'Geneid')
}


################################################################################
################################################################################
# Custom PCA plotting helper with some extra highlights

# Input vst/rlog normalized DESeq object. Created for instance with:
# DESeqDataSetFromMatrix(...) |> DESeq() |> vst()
x <- dat$`250` |> DESeq2::DESeq() |> DESeq2::vst()
ntop <- 500
color_values <- cbPalette[-5]
helper.pca <- function(x, plot.x = 'PC1', plot.y = 'PC2',
                       sh = 'photons', cl = 'CO2',
                       color_values = cbPalette[-5], ntop = 500) {
  # Adapted from DESeq2::plotPCA
  dba <- SummarizedExperiment::assay(x)
  rv <- MatrixGenerics::rowVars(dba)
  sel <- order(rv, decreasing = TRUE)[seq_len(min(ntop, length(rv)))]
  pca <- prcomp(t(dba[sel, ]))
  percentVar <- ( pca$sdev^2/sum(pca$sdev^2) ) %>%
    set_names(paste0('PC', 1:length(.)))
  
  x.meta <-
    x |>
    SummarizedExperiment::colData() |>
    as_tibble()
 
  pca$x |>
    as_tibble(rownames = 'lib') |>
    left_join(x.meta , 'lib') |>
    ggplot(aes(!! sym(x), !! sym(y),
               shape = !! sym(sh),
               color = !! sym(cl))) +
    geom_point(size = 5) +
    xlab(sprintf('%s: %.1f%% variance', plot.x, percentVar[[plot.x]] * 100)) +
    ylab(sprintf('%s: %.1f%% variance', plot.y, percentVar[[plot.y]] * 100)) +
    scale_color_manual(name = cl, values = color_values) +
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
