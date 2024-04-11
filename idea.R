#!/usr/bin/env Rscript

# idea on how to make better gene of interest selection
library(tidyverse)
library(clusterProfiler)

################################################################################

annot <-
  'data/C_annotation.tsv' |>
  read_tsv()

vst <-
  'analysis/D_vst-expression.tsv' |>
  read_tsv()

meta <-
  'data/C_meta.tsv' |>
  read_tsv() |>
  mutate_at('CO2', ~ fct_reorder(as.character(.x), as.numeric(.x)))

deg <-
  'analysis/D_stagewise-adjusted-DEGs.tsv' |>
  read_tsv()

raw.counts <- read_tsv('data/C_raw-counts.tsv')

raw.counts %>%
  select(- Geneid) %>%
  as.matrix() %>%
  magrittr::set_rownames(raw.counts$Geneid) -> raw.counts.mat

################################################################################

deg |>
  filter(is.de) |>
  select(Geneid) |>
  unique() |>
  left_join(vst) |>
  pivot_longer(- Geneid, names_to = 'lib') |>
  left_join(meta, 'lib') |>
  group_by(Geneid, CO2) |>
  summarize(avg = mean(value)) |>
  ggplot(aes(CO2, avg, group = Geneid)) +
  geom_line(alpha = .3)




################################################################################
################################################################################
################################################################################
# Sort of logFC
################################################################################
################################################################################
################################################################################
################################################################################
# Exclude rRNA

annot %>%
  filter(type != 'rRNA') %>%
  pull(Geneid) -> mask

raw.noribo.mat <- raw.counts.mat[mask, ]

################################################################################
# Run DESeq2+stageR

des <- DESeq2::DESeqDataSetFromMatrix(
  countData = raw.noribo.mat,
  colData = meta,
  design = ~ CO2
) |>
  DESeq2::DESeq()

# Contrast: 30% vs Rest
tribble(
  ~ coef,        ~ control, ~ case,
  "Intercept",       3,  1,
  "CO2_4_vs_0.04",   1,  0,
  "CO2_8_vs_0.04",   1,  0,
  "CO2_30_vs_0.04",  0,  1,
) |>
  mutate(
    # normalize by number of libs
    control  = control / 3,
    # create contrast
    contrast = case - control
  ) -> my.contrast

fc30 <-
  DESeq2::results(
    des,
    contrast = with(my.contrast, set_names(contrast, coef)),
    tidy = TRUE
  ) |>
  as_tibble() |>
  rename(Geneid = row)

################################################################################

fc30 |>
  with(plot(log2FoldChange, -log(padj)))

################################################################################
fc30 |>
  filter(padj <= 0.05) |>
  filter(abs(log2FoldChange) >= 1) |>
  left_join(annot) |>
  View()

fc30 |>
  filter(padj <= 0.05) |>
  filter(abs(log2FoldChange) >= 1) |>
  left_join(annot) |>
  drop_na(old_locus_tag) |>
  pull(old_locus_tag) |>
  write_lines('foo-string.txt')

################################################################################
lfc <-
  fc30 |>
  left_join(annot) |>
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
gsea |> as_tibble()


library(enrichplot)
i <- gsea
p <- gseaplot2(i, 1:nrow(i), 
               base_size = 14,
               color = RColorBrewer::brewer.pal(nrow(i), 'Paired'))

# Change axis text for clarity
p[[3]] <-  p[[3]] + ylab('log2 Fold Change')
# p[[1]] <- p[[1]] + ggtitle(j)

p

################################################################################
fisher.up <-
  fc30 |>
  filter(padj <= 0.05, log2FoldChange >= 1) |>
  left_join(annot) |>
  drop_na(old_locus_tag) |>
  pull(old_locus_tag) |>
  enrichMKEGG(
    organism = 'syp',
    maxGSSize = 70, # cutoff from comments above
    minGSSize = 5,
    pvalueCutoff = 0.05,
  ) |>
  as_tibble() |>
  mutate_at('Description', str_remove, ' - Picosynechococcus sp. PCC 7002')
fisher.up

################################################################################
fisher.down <-
  fc30 |>
  filter(padj <= 0.05, log2FoldChange <= -1) |>
  left_join(annot) |>
  drop_na(old_locus_tag) |>
  pull(old_locus_tag) |>
  enrichMKEGG(
    organism = 'syp',
    maxGSSize = 70, # cutoff from comments above
    minGSSize = 5,
    pvalueCutoff = 0.05,
  ) |>
  as_tibble() |>
  mutate_at('Description', str_remove, ' - Picosynechococcus sp. PCC 7002')

fisher.down
