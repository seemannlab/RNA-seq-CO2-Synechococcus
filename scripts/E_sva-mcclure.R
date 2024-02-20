#!/usr/bin/env Rscript
# Run DESeq after corection with sva

library(tidyverse)
library(ggpubr)
library(patchwork)

library(venn)
library(corrplot)

library(DESeq2)
library(sva)

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

meta.co2 <-
  'data/C_meta.tsv' |>
  read_tsv() |>
  mutate_at('batch', str_remove, '-Lane.$') |>
  mutate_at(c('CO2', 'photons'), function(x) {
    x |>
      as.character() |>
      fct_reorder(as.numeric(x))
  })

raw.co2 <-
  'data/C_raw-counts.tsv' |>
  read_tsv()

meta.mc <-
  'data/C_public-meta.tsv' |>
  read_tsv() |>
  filter(dataset == 'Photosynthetic-Growth-SRP063309') |>
  transmute(
    lib,
    CO2 = '0.04',
    photons = condition |>
      str_remove('mumol-') |>
      str_remove('photons') |>
      as.numeric(),
    batch = 'McClure et al.'
  ) |>
  mutate_at(c('CO2', 'photons'), function(x) {
    x |>
      as.character() |>
      fct_reorder(as.numeric(x))
  })

raw.mc <-
  'data/C_public-raw-counts.tsv' |>
  read_tsv() |>
  select(Geneid, !! meta.mc$lib)

annot <-
  'data/C_annotation.tsv' |>
  read_tsv()


################################################################################
# One combined dataset

meta <- meta.co2 |>
  mutate(batch = 'This study') |>
  bind_rows(meta.mc) |>
  mutate_at(c('CO2', 'photons'), function(x) {
    x |>
      as.character() |>
      fct_reorder(as.numeric(x))
  })

raw <- left_join(raw.co2, raw.mc, 'Geneid')

mat <-
  raw |>
  select(- Geneid) |>
  as.matrix() |>
  magrittr::set_rownames(raw$Geneid) |>
  # guarantee same order as rows in meta
  _[, meta$lib]

################################################################################
# Exclude rRNA and lowely expressed genes

annot %>%
  filter(type != 'rRNA') %>%
  pull(Geneid) -> mask

mat2 <- mat[mask, ]
mat2 <- mat2[rowMax(mat2) > 1, ]

################################################################################
# various versions of Combat correction

deseq.data <- list(
  'Without correction' = DESeqDataSetFromMatrix(
    countData = mat2,
    colData = meta,
    design = ~ CO2
  ),
  'With correction' = {
    x <- ComBat_seq(
      counts = mat2,
      # batch = with(meta, paste(batch, photons)),
      batch = meta$photons,
      group = NULL
    )
    DESeqDataSetFromMatrix(
      countData = x,
      colData = meta,
      design = ~ CO2
    )
  }
)

################################################################################
# Inspect impact of correction on normalization and PCA

deseq.des <- map(deseq.data, DESeq)

deseq.des |>
  map(sizeFactors) |>
  as_tibble() |>
  GGally::ggpairs() +
  theme_bw(18)

ggsave('analysis/E_sizeFactors.jpeg', width = 8, height = 8)
# -> DEseq2 normalization remains stable

################################################################################

deseq.vst <- map(deseq.data, vst)

my.plots <- function(des, path) {
  (
    helper.pca(des, sh = 'batch') /
      helper.pca(des, cl = 'photons',
             color_values = ggsci::pal_ucscgb()(8),
             sh = 'batch')
  ) |
    helper.skree(des)
  ggsave(path, width = 20, height = 10)
}

deseq.vst %>%
  map2(names(.), ~ my.plots(.x, sprintf(
    'analysis/E_pca-mcclure_%s.jpeg',
    .y |> str_replace_all('[^A-Za-z0-9]', '_')
  )))


################################################################################

des.hopes <-
  deseq.des |>
  map(function(xs) {
    # xs <- combats$`Photon AFTER dataset correction`
    des.hope <-
      DESeqDataSetFromMatrix(
        countData = counts(xs)[, meta.co2$lib],
        colData = meta.co2,
        design = ~ CO2
      )
    sizeFactors(des.hope) <- sizeFactors(xs)[meta.co2$lib]
    DESeq(des.hope)
  })



des.hopes |>
  map(vst) %>%
  map2(names(.), ~ my.plots(.x, sprintf(
    'analysis/E_pca-subset_%s.jpeg',
    .y |> str_replace_all('[^A-Za-z0-9]', '_')
  )))

################################################################################
# Impact on average logFC

tribble(
  # memo, low light was 1 and 15% condition, all others were high light
  ~ coef,        ~ low.light, ~ high.light,
  "Intercept",             2,  3 + 1,
  "CO2_1_vs_0.04",   1,  0,
  "CO2_4_vs_0.04",   0,  1,
  "CO2_8_vs_0.04",   0,  1,
  "CO2_15_vs_0.04",  1,  0,
  "CO2_30_vs_0.04",  0,  1,
) |>
  mutate(
    # normalize by number of libs
    low.light  = low.light / 2,
    high.light = high.light / 4,
    # create contrast
    contrast = high.light - low.light
  ) -> my.contrast

avs.light <-
  des.hopes |>
  map(~ results(
    .x,
    contrast = with(my.contrast, set_names(contrast, coef)),
    tidy = TRUE
  )) %>%
  map2(names(.), ~ mutate(.x, des = .y)) |>
  bind_rows()

################################################################################

avs.light |>
  filter(padj <= ALPHA) |>
  group_by(des) |>
  do(i = .$row) |>
  with(set_names(i, des)) |>
  venn::venn(zcolor = 'style', ilcs = 1.5, sncs = 1.5,
             box = FALSE, ggplot = TRUE) +
  annotate('text', 50, 950,
           label = 'logFC between light intensities, FDR 5%',
           hjust = 0, size = 8)
  
ggsave('analysis/E_venn-average-batch-effect.jpeg', width = 7, height = 7)


################################################################################

full.res <-
  des.hopes |>
  map(helper.pairwise.deg) %>%
  map2(names(.), function(res, i) {
    helper.stagewise(des.hopes[[i]], res, ALPHA) |>
      mutate(version = i)
  }) |>
  bind_rows()

write_tsv(full.res, 'analysis/E_stagewise-adjusted-DEGs.tsv')

################################################################################

full.res |>
  select(Geneid, test, log2FoldChange, version) |>
  spread(version, log2FoldChange) |>
  ggplot(aes(`With correction`, `Without correction`)) +
  geom_abline(slope = 1, color = 'blue') +
  geom_smooth(method = 'lm', color = 'red', se = FALSE) +
  geom_point(alpha = .7) +
  stat_cor(color = 'red') +
  facet_wrap(~ test, scales = 'free') +
  theme_bw(12)

ggsave('analysis/E_scatter-logFC-correction.jpeg', width = 12, height = 12)

################################################################################

full.res |>
  select(Geneid, test, log2FoldChange, version) |>
  spread(version, log2FoldChange) |>
  filter(abs(`With correction`) <= 5) |>
  ggplot(aes(`With correction`, `Without correction`)) +
  geom_abline(slope = 1, color = 'blue') +
  geom_smooth(method = 'lm', color = 'red', se = FALSE) +
  geom_point(alpha = .7) +
  stat_cor(color = 'red') +
  facet_wrap(~ test, scales = 'free') +
  theme_bw(12)

ggsave('analysis/E_scatter-logFC-correction-outliers.jpeg', width = 12, height = 12)

################################################################################

full.res |>
  filter(padj.stageR <= ALPHA) |>
  select(version, Geneid) |>
  unique() |>
  group_by(version) |>
  do(i = .$Geneid) |>
  with(set_names(i, version)) |>
  venn::venn(zcolor = 'style', ilcs = 1.5, sncs = 1.5,
             box = FALSE, ggplot = TRUE) +
  annotate('text', 50, 950,
           label = 'stage-wise adjusted FDR 5%',
           hjust = 0, size = 8)

ggsave('analysis/E_venn-deg.jpeg', width = 7, height = 7)

################################################################################

des.hopes |>
  names() |>
  map(function(xs) {
    # xs <- 'With correction'
    full.res |>
      filter(version == xs, is.de) |>
      pull(Geneid) |>
      unique() -> mask
    
    cl <- list(
      'CO2' = meta.co2 %>%
        pull(CO2) |>
        levels() %>%
        set_names(cbPalette[-5][1:(length(.))], .)
    )
    
    with(
      meta.co2,
      data.frame('CO2' = as.character(CO2), row.names = lib)
    ) -> cl.df
    
    norm.vst <-
      des.hopes[[xs]] |>
      DESeq2::vst() |>
      SummarizedExperiment::assay()
    # Cluster columns ahead of heatmap to rotate tree nicely
    dat <-  norm.vst[mask, ]
    
    lib.clust <-
      dat |>
      t() |>
      dist() |>
      hclust() |>
      ape::as.phylo() |>
      ape::rotateConstr(
        meta.co2 |>
          arrange(CO2, sample) |>
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
      filename = sprintf(
        'analysis/E_-heatmap-%s.jpeg',
        xs |> str_replace_all('[^A-Za-z0-9]', '_')
      )
    )
  })
dev.off()


