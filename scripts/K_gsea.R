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

path.gene <-
  'analysis/I_gene2pathway.tsv' |>
   read_tsv()

vst <-
  'analysis/D_vst-expression.tsv' |>
  read_tsv()

################################################################################
# Only check the logFC to next higher CO2 concentration

# todo <-
#   meta |>
#   pull(CO2) |>
#   unique() |>
#   sort() %>%
#   tibble(i = ., j = lead(.)) |>
#   drop_na() |>
#   with(sprintf('%s->%s%% CO2', i, j))

# test for all
todo <-
  deg |>
  pull(test) |>
  unique()

################################################################################
# Determine max size cutoff to exclude overly large pathways
#
# path.gene |>
#   count(Pathway, Title) |>
#   arrange(desc(n))
# 
# syp01100 Metabolic pathways                             604   \
# syp01110 Biosynthesis of secondary metabolites          264    |
# syp01240 Biosynthesis of cofactors                      127    | exclude
# syp01120 Microbial metabolism in diverse environments   117    |
# syp01230 Biosynthesis of amino acids                     89    |
# syp02010 ABC transporters                                76   /
# syp00970 Aminoacyl-tRNA biosynthesis                     67
# syp00195 Photosynthesis                                  66
# syp01200 Carbon metabolism                               66
# syp03010 Ribosome                                        58

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
      maxGSSize = 70, # cutoff from comments above
      minGSSize = 5,
      # pvalueCutoff = 0.01,
      pvalueCutoff = 0.05,
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
# Overview all enrichment results from all test in a single figure

my.order <- c(
  "0.04->4% CO2",
  "0.04->8% CO2",
     "4->8% CO2",
     "4->30% CO2", 
     "8->30% CO2",
  "0.04->30% CO2"
) |>
  str_remove(' CO2') |>
  str_replace('->', '\U2192')

res |>
  mutate_at('test', str_remove, ' CO2') |>
  mutate_at('test', str_replace, '->', '\U2192') |>
  mutate_at('test', fct_relevel, my.order) |>
  ggplot(aes(test, fct_rev(Description), fill = NES)) +
  geom_tile() +
  scale_fill_gradient2(high = 'red', low = 'blue',
                       name = 'Enrichment Score') +
  ylab(NULL) +
  xlab(NULL) +
  theme_pubr(18) +
  theme(axis.text.x = element_text(angle = 60, hjust = 1))

ggsave('analysis/K_enrichment-overview.jpeg',
       width = 10, height = 10, dpi = 400)

################################################################################
# z-scale expression

vst.mat <-
  vst |>
  select(- Geneid) |>
  as.matrix() |>
  magrittr::set_rownames(vst$Geneid)

z.mat <-
  vst.mat |>
  apply(1, scale) |>
  t() |> 
  magrittr::set_colnames(colnames(vst.mat))

z.expr <-
  z.mat |>
  as_tibble(rownames = 'Geneid') |>
  pivot_longer(- Geneid, names_to = 'lib', values_to = 'z.expression') |>
  left_join(meta, 'lib')

################################################################################
# Overview of expression patterns



path.dat <-
  res |>
  select(Pathway = ID, Description) |>
  unique() |>
  left_join(path.gene, 'Pathway') |>
  select(Pathway = Description, Geneid) |>
  semi_join(
    deg |>
      filter(is.de) |>
      select(Geneid) |>
      unique(),
    'Geneid'
  ) |>
  left_join(z.expr, 'Geneid', relationship = "many-to-many")

# Pretty plot please
path.dat |>
  mutate_at('Pathway', str_replace, '(?<=.{15}) ', '\n') |>
  mutate_at('CO2', ~ fct_reorder(as.character(.x), .x)) |>
  ggplot(aes(CO2, z.expression)) +
  geom_violin(aes(fill = Pathway), show.legend = FALSE) +
  # scale_fill_brewer(palette = 'Spectral') +
  ggsci::scale_fill_ucscgb() +
  geom_boxplot(width = .2) +
  geom_smooth(aes(x = as.integer(CO2), group = 1),
              method = 'loess', se = FALSE,
              color = 'red') +
  facet_wrap(~ Pathway, ncol = 4) +
  ylab('z-scaled log expression') +
  xlab(NULL) +
  theme_pubr(18)

ggsave('analysis/K_expression-overview.jpeg',
       width = 12, height = 12)

################################################################################


sessionInfo()
