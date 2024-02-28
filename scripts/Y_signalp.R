#!/usr/bin/env Rscript

# Inspect the signalP results

library(tidyverse)
library(ggpubr)
library(patchwork)

library(conflicted)
conflicts_prefer(dplyr::filter)
conflicts_prefer(dplyr::rename)
conflicts_prefer(purrr::reduce)
conflicts_prefer(dplyr::count)

# https://riptutorial.com/r/example/28354/colorblind-friendly-palettes
cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442",
               "#0072B2", "#D55E00", "#CC79A7")

################################################################################
################################################################################
# Load data

annot <-
  'data/C_annotation.tsv' |>
  read_tsv()

deg <-
  'analysis/E_dge-stagewise-analysis.tsv' |>
  read_tsv()

meta <-
  'data/C_meta.tsv' |>
  read_tsv()

vst <-
  'analysis/E_vst.tsv' |>
  read_tsv()

signal <-
  'analysis/J_signalp/prediction_results.txt' |>
  read_tsv(skip = 1) |>
  rename('Geneid' = '# ID')

################################################################################
# nicer names for the signals
signal.names <- c(
  'OTHER' = 'No SP',
  'SP' = 'Sec/SPI',
  'LIPO' = 'Sec/SPII',
  'TAT' = 'Tat/SPI',
  'TATLIPO' = 'Tat/SPII',
  'PILIN' = 'Sec/SPIII'
)

signal |>
  mutate(Prediction = signal.names[Prediction]) |>
  count(Prediction)
# Prediction     n
# No SP       2935
# Sec/SPI      166
# Sec/SPII      67
# Sec/SPIII      4
# Tat/SPI        9
# Tat/SPII       5

################################################################################

signal |>
  select(- Prediction, - `CS Position`) |>
  pivot_longer(- Geneid) |>
  # group_by(Geneid) |> summarize(x = sum(value)) |> pull(x) |> summary()
  # values per gene add up to 1 ± -.0001
  mutate(
    name = str_remove(name, '\\(.*$'),
    name = signal.names[name]
  ) |>
  ggplot(aes(value, color = name)) +
  stat_ecdf(size = 1.2) +
  # ggsci::scale_color_jama(name = 'Signal peptide') +
  scale_color_manual(
    name = 'SingalP prediction',
    values = signal |>
      select(Prediction) |>
      unique() |>
      nrow() |>
      RColorBrewer::brewer.pal('Paired')
  ) +
  scale_x_continuous(breaks = seq(0, 1, .1)) +
  scale_y_continuous(breaks = seq(0, 1, .1)) +
  xlab('Signal peptide probability') +
  ylab('Cumulative density function') +
  theme_pubr(18) 

ggsave('analysis/J_signal-probs.jpeg',
       width = 10, height = 7, dpi = 400)

################################################################################
# Heatmap of expression for genes with signals

vst.mat <-
  vst |>
  select(- Geneid) |>
  as.matrix() |>
  magrittr::set_rownames(vst$Geneid)

row.annot <-
  signal |>
  filter(Prediction != 'OTHER') |>
  mutate(Prediction = signal.names[Prediction]) |>
  with(data.frame('SignalP' = Prediction,
                  check.names = FALSE,
                  row.names = Geneid))

col.annot <-
  meta |>
  with(data.frame('CO2' = as.character(condition), row.names = lib))

annot.cl <- list(
  'SignalP' = row.annot |>
    select(SignalP) |>
    unique() |>
    arrange(SignalP) %>%
    # mutate(cl = ggsci::pal_jama()(nrow(.) + 1)[-1]) |>
    mutate(cl = RColorBrewer::brewer.pal(nrow(.) + 1, 'Paired')[- 1]) |>
    with(set_names(cl, SignalP)),
  'CO2' =  meta |>
    select(CO2 = condition) |>
    unique() |>
    arrange(CO2) |>
    mutate(cl = cbPalette[-5][1:n()]) |>
    with(set_names(cl, as.character(CO2)))
)

# only on the genes with signal
dat <- vst.mat[rownames(row.annot), ]

# Cluster columns ahead of pheatmap to rotate tree nicely
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
  show_rownames = FALSE,
  show_colnames = FALSE,
  cluster_cols = lib.clust,
  annotation_row = row.annot,
  annotation_col = col.annot,
  annotation_colors = annot.cl,
  fontsize = 12,
  color = colorRampPalette(rev(
    RColorBrewer::brewer.pal(n = 7, name = "RdBu")))(59),
  filename = 'analysis/J_heatmap.jpeg'
)
dev.off()

################################################################################

signal |>
  filter(Prediction != 'OTHER') |>
  select(Geneid, Prediction) |>
  left_join(annot, 'Geneid') -> foo

bar <-
  deg |>
  filter(is.de, test %in% c('4->30% CO2', '8->30% CO2')) |>
  select(Geneid, test, log2FoldChange) |>
  mutate(test = paste('logFC', test)) |>
  spread(test, log2FoldChange)

baz <-
  foo |>
  filter(! str_detect(product, 'hypothetical')) |>
  filter(! str_detect(product, 'domain-containing protein')) |>
  # left_join(bar, 'Geneid') |>
  inner_join(bar, 'Geneid') |>
  select(-Geneid) |>
  mutate_at('Prediction', ~ signal.names[.x]) |>
  select(
    contains('locus'),
    name, product,
    signalP.prediction = Prediction,
    contains('logFC')
  ) |>
  arrange(`logFC 4->30% CO2`)

write_tsv(baz, 'analysis/J_signalp-subset.tsv')

################################################################################
# Enrichment of signal in pathways?

# Collect loci per signals
loci.sets <-
  signal |>
  filter(Prediction != 'OTHER') |>
  select(Geneid, Prediction) |>
  mutate(Prediction = signal.names[Prediction]) |>
  left_join(annot, 'Geneid') |>
  drop_na(old_locus_tag) |>
  select(Prediction, old_locus_tag) |>
  unique() |>
  group_by(Prediction) |>
  do(i = list(.$old_locus_tag)) |>
  with(set_names(i, Prediction)) |>
  map(1)

# Enrichment test
signal.enrich <-
  loci.sets |>
  map(
    clusterProfiler::enrichKEGG,
    organism = 'syp',
    pvalueCutoff = 0.05,
    pAdjustMethod = "fdr"
  ) |>
  # combine tables
  map(as_tibble) %>%
  map2(names(.), ~ mutate(.x, signal = .y)) |>
  bind_rows() %>%
  transmute(
    signal, 
    KEGG = ID,
    Pathway = str_remove(Description, ' - Picosynechococcus sp. PCC 7002$'),
    GeneRatio, BgRatio,
    Pvalue = pvalue,
    FDR = p.adjust,
    genes = geneID
  ) |>
  # Extract enrichment factor
  mutate(
    path.size = BgRatio |>
      str_remove('/.*') |>
      as.integer(),
    genes.with.signal.in.path = GeneRatio |>
      str_remove('/.*') |>
      as.integer(),
    genes.with.signal = GeneRatio |>
      str_remove('.*/') |>
      as.integer(),
    prop.all.genes.in.pathway = BgRatio |>
      map(~ parse(text = .x)) |>
      map(eval) |>
      unlist(),
    expected.genes = genes.with.signal * prop.all.genes.in.pathway,
    enrichment = genes.with.signal.in.path / expected.genes
  )

write_tsv(signal.enrich, 'analysis/J_signal-enrichment.tsv')

################################################################################

signal.enrich |>
  mutate(
    y = sprintf('%s - %s (%s)', signal, Pathway, KEGG) |>
      fct_reorder(enrichment),
    nice = sprintf(
      'Pathway size: %g\nGenes with signal: %g\nEnrichment: %.1f, FDR: %.1e',
      path.size,
      genes.with.signal.in.path, 
      enrichment, FDR
    ),
    hj = ifelse(enrichment > 20, 1, 0),
    x2 = enrichment + ifelse(enrichment > 20, -1, 1) * 1.5,
  ) |>
  ggplot(aes(enrichment, y, color = - log10(FDR),
             size = path.size, label = nice)) +
  # geom_point(size = 8) +
  geom_point() +
  scale_size(range = c(5, 20), name = 'Pathway size') +
  geom_text(aes(x = x2, hjust = hj), color = 'black', size = 5) +
  scale_color_viridis_c() +
  guides(color = guide_colorbar(barwidth = unit(5, 'cm'))) +
  ylab(NULL) +
  xlab('Enrichment ratio observed over expected genes in pathway') +
  theme_pubr(18) +
  theme(
    panel.grid.major.y = element_line(linetype = 'dotted')
  )

ggsave('analysis/J_signal-enrichment.jpeg',
       width = 14, height = 8, dpi = 400)

################################################################################
# Query pathway associations

KEGG_DATA <- clusterProfiler:::prepare_KEGG('syp', 'KEGG', 'kegg')
tibble(
  Pathway = names(KEGG_DATA$PATHID2NAME),
  Title = str_remove(KEGG_DATA$PATHID2NAME, ' - Picosynechococcus sp. PCC 7002$')
) %>%
  left_join(
    KEGG_DATA$PATHID2EXTID %>%
      map(~ tibble(old_locus_tag = .x)) %>%
      map2(names(.), ~ mutate(.x, Pathway = .y)) %>%
      bind_rows(),
    multiple = "all"
  ) %>%
  left_join(annot, 'old_locus_tag') -> gene.path


write_tsv(gene.path, 'analysis/J_gene2pathway.tsv')

################################################################################

signal.enrich |>
  select(signal, KEGG, Pathway) |>
  left_join(gene.path, c('KEGG' = 'Pathway', 'Pathway' = 'Title'),
            relationship = "many-to-many") |>
  left_join(signal, 'Geneid') |>
  filter(signal == signal.names[Prediction]) |>
  select(signal, KEGG, Pathway, old_locus_tag, Geneid, name, product) |>
  write_tsv('analysis/J_enriched-genes.tsv')

