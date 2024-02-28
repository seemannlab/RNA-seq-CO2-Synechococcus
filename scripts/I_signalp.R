#!/usr/bin/env Rscript

# Inspect the signalP results

library(tidyverse)
library(ggpubr)
library(patchwork)

library(clusterProfiler)

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
  'analysis/D_stagewise-adjusted-DEGs.tsv' |>
  read_tsv()

meta <-
  'data/C_meta.tsv' |>
  read_tsv()

vst <-
  'analysis/D_vst-expression.tsv' |>
  read_tsv()

signal <-
  'analysis/H_signalp/prediction_results.txt' |>
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

################################################################################
# density of prediciton scores

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

ggsave('analysis/I_signal-probs.jpeg',
       width = 10, height = 7, dpi = 400)

################################################################################
# make prettier table

signal2 <-
  signal |>
  transmute(
    Geneid,
    signal = signal.names[Prediction]
  ) |>
  left_join(annot, 'Geneid') |>
  left_join(
    deg |>
      filter(is.de) |>
      transmute(Geneid, is.de) |>
      unique(),
    'Geneid'
  ) |>
  mutate_at('is.de', replace_na, FALSE)

write_tsv(signal2, 'analysis/I_signals.tsv')

signal2 |>
  count(signal, is.de) |>
  spread(is.de, n, fill = 0) |>
  transmute(
    signal,
    genes.total = `TRUE` + `FALSE`,
    diff.expressed = `TRUE`,
    ratio = diff.expressed / genes.total
  ) |>
  write_tsv('analysis/I_overview.tsv')


################################################################################
# Heatmap of expression for genes with signals

vst.mat <-
  vst |>
  select(- Geneid) |>
  as.matrix() |>
  magrittr::set_rownames(vst$Geneid)

row.annot <-
  signal2 |>
  filter(signal != 'No SP', is.de) |>
  with(data.frame('SignalP' = signal,
                  check.names = FALSE,
                  row.names = Geneid))

col.annot <-
  meta |>
  with(data.frame('CO2' = as.character(CO2), row.names = lib))

annot.cl <- list(
  'SignalP' = row.annot |>
    select(SignalP) |>
    unique() |>
    arrange(SignalP) %>%
    # mutate(cl = ggsci::pal_jama()(nrow(.) + 1)[-1]) |>
    mutate(cl = RColorBrewer::brewer.pal(nrow(.) + 1, 'Paired')[- 1]) |>
    with(set_names(cl, SignalP)),
  'CO2' =  meta |>
    select(CO2) |>
    unique() |>
    arrange(CO2) |>
    mutate(cl = cbPalette[c(1, 6, 2, 7)]) |>
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
      arrange(CO2, sample) |>
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
  filename = 'analysis/I_heatmap.jpeg'
)
dev.off()

################################################################################
# Enrichment of signal in pathways?

# Collect loci per signals
loci.sets <-
  signal2 |>
  filter(signal != 'No SP', is.de) |>
  select(Geneid, signal) |>
  left_join(annot, 'Geneid') |>
  drop_na(old_locus_tag) |>
  select(signal, old_locus_tag) |>
  unique() |>
  group_by(signal) |>
  do(i = .$old_locus_tag) |>
  with(set_names(i, signal))

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

write_tsv(signal.enrich, 'analysis/I_signal-enrichment.tsv')

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

ggsave('analysis/I_signal-enrichment.jpeg',
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


write_tsv(gene.path, 'analysis/I_gene2pathway.tsv')

################################################################################

signal.enrich |>
  select(signal, KEGG, Pathway) |>
  left_join(gene.path, c('KEGG' = 'Pathway', 'Pathway' = 'Title'),
            relationship = "many-to-many") |>
  left_join(signal, 'Geneid') |>
  filter(signal == signal.names[Prediction]) |>
  select(signal, KEGG, Pathway, old_locus_tag, Geneid, name, product) |>
  write_tsv('analysis/I_enriched-genes.tsv')


################################################################################
sessionInfo()