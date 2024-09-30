#!/usr/bin/env Rscript

# Purpose: GSEA of all pairwise logFCs

library(tidyverse)

library(cowplot)
library(patchwork)
library(ggpubr)

library(clusterProfiler)
library(enrichplot)

# https://riptutorial.com/r/example/28354/colorblind-friendly-palettes
cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442",
               "#0072B2", "#D55E00", "#CC79A7")

################################################################################
# Load input data

annot <-
  'data/C_annotation.tsv' |>
  read_tsv()

deg <-
  'analysis/D_stagewise-adjusted-DEGs.tsv' |>
  read_tsv()


################################################################################
# Compute geneset enrichment of KEGG pathways per logFC


helper <- function(x) {
  # x <- '0.04->30% CO2'
  gsea <-
    deg |>
    filter(test == x) |>
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
  
  
  i <- gsea
  p <- gseaplot2(i, 1:nrow(i), 
                 base_size = 14,
                 color = RColorBrewer::brewer.pal(nrow(i), 'Paired'))
  
  # Change axis text for clarity
  p[[3]] <-  p[[3]] + ylab('log2 Fold Change')
  
  gsea |>
    as_tibble() |>
    write_tsv(sprintf('analysis/M_%s.tsv', x |> str_replace_all('[^a-zA-Z0-9]', '_')))
  
  ggsave(sprintf('analysis/M_%s.jpeg', x |>
                   str_replace_all('[^a-zA-Z0-9]', '_')),
         plot = p,
         width = 12, height = 12, dpi = 400)
}

################################################################################
deg |> pull(test) |> unique() |> map(helper)

################################################################################
