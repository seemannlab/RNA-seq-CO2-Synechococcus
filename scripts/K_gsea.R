# Inspect GSEA per logFC for potential pathways
# (Alternative to fisher test from stringApp)

set.seed(42)

library(tidyverse)
library(ggpubr)
library(patchwork)

library(clusterProfiler)
library(enrichplot)

library(conflicted)
conflicts_prefer(dplyr::filter)
conflicts_prefer(dplyr::rename)
conflicts_prefer(purrr::reduce)
conflicts_prefer(dplyr::count)

# https://riptutorial.com/r/example/28354/colorblind-friendly-palettes
cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442",
               "#0072B2", "#D55E00", "#CC79A7")

################################################################################
# Load data

annot <-
  'data/C_annotation.tsv' |>
  read_tsv()

deg <-
  'analysis/D_stagewise-adjusted-DEGs.tsv' |>
  read_tsv()

meta <-
  'data/C_meta.tsv' |>
  read_tsv()

################################################################################
# Only check the logFC to next higher CO2 concentration

todo <-
  meta |>
  pull(CO2) |>
  unique() |>
  sort() %>%
  tibble(i = ., j = lead(.)) |>
  drop_na() |>
  with(sprintf('%s->%s%% CO2', i, j))

################################################################################
# Compute GSEA

gs.res <-
  todo |>
  set_names(todo) |>
  map(function(i) {
    # i <- "8->30% CO2"
    # Vector of logFCs
    lfc <-
      deg |>
      filter(test == i) |>
      drop_na(old_locus_tag, log2FoldChange) |>
      arrange(desc(log2FoldChange)) |>
      with(set_names(log2FoldChange, old_locus_tag))
    
    # Compute geneset enrichment of KEGG pathways
    gsea <- gseKEGG(
      geneList = lfc,
      organism = 'syp',
      pvalueCutoff = 0.01,
      eps = 0
    )
    gsea@result$Description <- str_remove(gsea$Description, ' - Picosynechococcus sp. PCC 7002')
    gsea
  })

################################################################################
# Pretty plots

gs.res %>%
  map2(names(.), function(i, j) {
    p <- gseaplot2(i, 1:nrow(i), 
                   base_size = 14,
                   color = RColorBrewer::brewer.pal(nrow(i), 'Paired'))
    
    # Change axis text for clarity
    p[[3]] <-  p[[3]] + ylab('log2 Fold Change')
    p[[1]] <- p[[1]] + ggtitle(j)
    ggsave(sprintf(
      'analysis/K_gsea-%s.jpeg',
      str_replace_all(j, '[^a-zA-Z0-9-_]', '.')
      ), 
      plot = p,
      width = 9, height = 12, dpi = 400
    )
  })

################################################################################

res <-
  gs.res |>
  map(as_tibble) |>
  map(select, - leading_edge, - core_enrichment) %>%
  map2(names(.), ~ mutate(.x, test = .y)) |>
  bind_rows()

write_tsv(res, 'analysis/K_gsea.tsv')


################################################################################
sessionInfo()
