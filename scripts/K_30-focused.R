#!/usr/bin/env Rscript

# Purpose: Expression analysis, but focused on differences to 30%

library(tidyverse)
library(gt)

library(cowplot)
library(patchwork)
library(ggpubr)

library(enrichplot)
library(clusterProfiler)

# For saving GT table to temporary png
Sys.setenv(CHROMOTE_CHROME = '/Applications/Brave Browser.app/Contents/MacOS/Brave Browser')

# https://riptutorial.com/r/example/28354/colorblind-friendly-palettes
cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442",
               "#0072B2", "#D55E00", "#CC79A7")

################################################################################
# Load input data

annot <-
  'data/C_annotation.tsv' |>
  read_tsv()

meta <-
  'data/C_meta.tsv' |>
  read_tsv() |>
  mutate_at('CO2', ~ fct_reorder(as.character(.x), as.numeric(.x)))

raw.counts <-
  'data/C_raw-counts.tsv' |>
  read_tsv()

deg <-
  'analysis/D_stagewise-adjusted-DEGs.tsv' |>
  read_tsv()

vst <-
  'analysis/D_vst-expression.tsv' |>
  read_tsv()


################################################################################
# Build matrices

vst.mat <-
  vst |>
  select(- Geneid) %>%
  as.matrix() %>%
  magrittr::set_rownames(vst$Geneid)

raw.counts.mat <-
  raw.counts %>%
  select(- Geneid) %>%
  as.matrix() %>%
  magrittr::set_rownames(raw.counts$Geneid)

################################################################################
# Exclude rRNA

annot %>%
  filter(type != 'rRNA') %>%
  pull(Geneid) -> mask

raw.noribo.mat <- raw.counts.mat[mask, ]

################################################################################
# Run DESeq2 for 30% against all other conditions

meta2 <-
  meta |>
  mutate(CO2 = ifelse(CO2 == 30, '30', 'other') |>
           fct_relevel('other', '30'))

des <- DESeq2::DESeqDataSetFromMatrix(
  countData = raw.noribo.mat,
  colData = meta2,
  design = ~ CO2
) |>
  DESeq2::DESeq()


deg30 <-
  DESeq2::results(
    des,
    contrast = c('CO2', '30', 'other'),
    tidy = TRUE
  ) |>
  as_tibble() |>
  rename(Geneid = row)

write_tsv(deg30, 'analysis/K_logFC-vs-30.tsv')

################################################################################
deg |>
  select(Geneid, test, log2FoldChange.pairwise = log2FoldChange) |>
  left_join(deg30, 'Geneid') |>
  ggscatter(
    'log2FoldChange', 'log2FoldChange.pairwise',
    add = 'reg.line', add.params = list(color = 'red'),
    cor.coef = TRUE, cor.coeff.args = list(color = 'red', size = 5),
    facet.by = 'test',
    alpha = .5
  ) +
  xlab('30% CO2 vs all other conditions') +
  ylab('Pairwise comparison') +
  theme_pubr(18) +
  ggtitle('log2 Fold-Changes')

ggsave('analysis/K_logFC-comparisons.jpeg',
       width = 8, height = 8, dpi = 400)

################################################################################
# Table for logFC / P-Value cutoff

deg30 |>
  mutate_at('padj', cut,
            c(.1, .05, .01, .001, 0),
            include.lowest = TRUE) |>
  mutate_at('log2FoldChange', cut,
            c(-Inf, -2, -1, -.5, .5, 1, 2, Inf)) |>
  count(log2FoldChange, padj) |>
  drop_na() |>
  spread(padj, n, fill = 0) |>
  gt(rowname_col = 'log2FoldChange') |>
  # fmt_bins(sep = ' - ', ) |>
  text_transform(
    \(x) str_replace(x, 'Inf', '∞'),
    locations = cells_stub()
  ) |>
  tab_spanner(
    label = "FDR adjusted P-Value",
    columns = everything()
  ) |>
  tab_stubhead(label = 'log2 Fold-Change') |>
  data_color(palette = "viridis") |>
  gtsave('tmp-foo.png', zoom = 5)


# load as grid
p.tab <-
  ggdraw() +
  draw_image('tmp-foo.png') +
  theme_pubr(18) +
  theme(axis.line = element_blank())

file.remove('tmp-foo.png')

################################################################################
# Heatmap of expression

deg30 %>%
  filter(padj <= 0.001) |>
  filter(abs(log2FoldChange) >= 1) |>
  pull(Geneid) %>%
  unique -> mask

cl <- list(
  'CO2' = meta %>%
    pull(CO2) %>%
    levels %>%
    set_names(cbPalette[c(1, 6, 2, 7)], .)
)

with(
  meta,
  data.frame('CO2' = as.character(CO2), row.names = lib)
) -> cl.df

# Cluster columns ahead of pheatmap to rotate tree nicely
dat <-  vst.mat[mask, ]
# z-scale data ahead of heatmap to keep clustering of rows consistent
dat <-
  dat |>
  apply(1, scale) |>
  t() |>
  magrittr::set_colnames(colnames(dat))

lib.clust <-
  dat |>
  t() |>
  dist() |>
  hclust() |>
  ape::as.phylo() |>
  ape::rotateConstr(
    meta |>
      arrange(CO2, sample) |>
      pull(lib)
  ) |>
  as.hclust()

pheatmap::pheatmap(
  dat,
  border_color = NA,
  scale = 'none',
  cluster_cols = lib.clust,
  show_rownames = FALSE,
  show_colnames = FALSE,
  annotation_col = cl.df,
  annotation_colors = cl,
  color = colorRampPalette(rev(
    RColorBrewer::brewer.pal(n = 7, name = "RdBu")))(59),
  filename = 'analysis/K_heatmap-30-focused.jpeg'
)
dev.off()

# load as grid
p.heat <-
  ggdraw() +
  draw_image(
    'analysis/K_heatmap-30-focused.jpeg'
  ) +
  theme_pubr(18) +
  theme(axis.line = element_blank())



################################################################################
# Compute geneset enrichment of KEGG pathways

gsea <-
  deg30 |>
  left_join(annot) |>
  drop_na(old_locus_tag, log2FoldChange) |>
  arrange(desc(log2FoldChange)) |>
  with(set_names(log2FoldChange, old_locus_tag)) |>
  gseKEGG(
    organism = 'syp',
    maxGSSize = 70, # cutoff from comments above
    minGSSize = 5,
    pvalueCutoff = 0.05,
    eps = 0
  )
gsea@result$Description <- str_remove(gsea$Description, ' - Picosynechococcus sp. PCC 7002')

gsea |>
  as_tibble() |>
  write_tsv('analysis/K_gsea.tsv')


i <- gsea
p <- gseaplot2(i, 1:nrow(i), 
               base_size = 14,
               color = RColorBrewer::brewer.pal(nrow(i), 'Paired'))

# Change axis text for clarity
p[[3]] <-  p[[3]] + ylab('log2 Fold Change')

################################################################################

p2 <-
  wrap_elements(full = print(p)) +
 theme_pubr(18)

(p.tab | p.heat) / p2 +
  plot_layout(widths = c(1, 1.5), heights = c(1, 2)) +
  plot_annotation(tag_levels = 'A', title = '30% vs all other conditions')

ggsave('analysis/K_overview-30-focused.jpeg',
       width = 12, height = 12, dpi = 400)

################################################################################

annot |>
  filter(Geneid %in% mask) |>
  drop_na(old_locus_tag) |>
  pull(old_locus_tag) |>
  write_lines('analysis/K_string-30-focused.txt')
