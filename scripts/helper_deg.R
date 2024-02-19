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
    ggplot(aes(!! sym(plot.x), !! sym(plot.y),
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
################################################################################

helper.skree <- function(x, ntop = 500) {
  # Adapted from DESeq2::plotPCA
  dba <- SummarizedExperiment::assay(x)
  rv <- MatrixGenerics::rowVars(dba)
  sel <- order(rv, decreasing = TRUE)[seq_len(min(ntop, length(rv)))]
  pca <- prcomp(t(dba[sel, ]))
  percentVar <- ( pca$sdev^2/sum(pca$sdev^2) ) %>%
    set_names(paste0('PC', 1:length(.)))
  
  dat <-
    percentVar %>%
    tibble(
      pc = names(.) |> str_remove('PC') |> as.integer(),
      'Explained varinace' = . * 100,
      'Cumulative' = cumsum(.) * 100
    ) |>
    select(- `.`) |>
    pivot_longer(- pc) |>
    mutate(nice = ifelse(
      (pc <= 4) & (pc != 1 | name == 'Cumulative'),
      sprintf('%.1f', value),
      NA_character_
    ))
  dat |>
    ggplot(aes(pc, value, label = nice,
               group = name, color = name)) +
    geom_line() +
    geom_point(size = 3) +
    xlab('Principle component') +
    ylab('Variance %') +
    ggsci::scale_color_jama(name = NULL) +
    # change repel y nudge for curves
    ggrepel::geom_text_repel(
      nudge_x = 5, nudge_y = -4,
      show.legend = FALSE,
      size = 5, color = 'black',
      data = filter(dat, name == 'Cumulative')
    ) +
    ggrepel::geom_text_repel(
      nudge_x = 5, nudge_y = +4,
      show.legend = FALSE,
      size = 5, color = 'black',
      data = filter(dat, name != 'Cumulative')
    ) +
    theme_pubr(18)
}

################################################################################
