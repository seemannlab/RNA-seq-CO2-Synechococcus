#!/usr/bin/env Rscript

## NOTE: Cytoscape with the network session loaded needs to be running for this script

# Overview of the StringApp generated cluster in relationship to the
# logFC from the DEG analysis


library(tidyverse)
library(ggpubr)
library(patchwork)

library(mixOmics)

library(conflicted)

conflicts_prefer(purrr::map)
conflicts_prefer(dplyr::select)
conflicts_prefer(dplyr::filter)
conflicts_prefer(dplyr::count)

# https://riptutorial.com/r/example/28354/colorblind-friendly-palettes
cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442",
               "#0072B2", "#D55E00", "#CC79A7")

library(furrr)
plan(multisession)

################################################################################

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

vst.mat <-
  vst |>
  select(- Geneid) |>
  as.matrix() |>
  magrittr::set_rownames(vst$Geneid)

################################################################################
# Query clustering info per genes

node.dat <- RCy3::getTableColumns()

mcl <-
  node.dat |>
  with(tibble(locus = name, string.name = `display name`, cluster = `__mclCluster`)) |>
  drop_na(cluster) |>
  mutate(locus = str_remove(locus, '32049.')) |>
  unique() |>
  left_join(annot, c('locus' = 'old_locus_tag'))

mcl |>
  write_tsv('analysis/L_mcl-clustering.tsv')
# mcl <- read_tsv('analysis/L_mcl-clustering.tsv')

################################################################################
# use command file to export all the clusters enrichment from cytoscape

x <- 1:100

# Query cytoscape for enrichment per cluster
c(
  'string retrieve group-wise enrichment maxGroupNumber=100 minGroupSize=5',
  sprintf(
    'table export outputFile="%s/analysis/L_cluster-enrichment/cluster_%g.csv" table="STRING Enrichment: __mclCluster %g"',
    getwd(), x, x
  )
) |>
  write_lines('analysis/L_cluster-enrichment-export-commands.txt')

dir.create('analysis/L_cluster-enrichment/')
tryCatch({
  RCy3::commandRunFile(file.path(
    getwd(),
    'analysis/L_cluster-enrichment-export-commands.txt'
  ))
},
# ignore error on "empyt tables'
error=function(error_message) {}
)

# load the exported data
string.enrichment <-
  'analysis/L_cluster-enrichment/*csv' |>
  Sys.glob() %>%
  set_names(. |> basename() |> str_remove('.csv$')) |>
  map(read_csv) %>%
  map2(names(.), ~ mutate(.x, cluster = .y)) |>
  bind_rows()

# Cleanup string enrichment export
enrich <-
  string.enrichment |>
  select(
    cluster,
    database = category,
    geneset = `term name`,
    geneset.name = description,
    no.genes = `# genes`,
    no.background = `# background genes`,
    pvalue = `p-value`,
    FDR = `FDR value`
  ) |>
  left_join(
    mcl |>
      mutate(cluster = paste0('cluster_', cluster)) |>
      count(cluster, name = 'cluster.size'),
    'cluster'
  )

# match display name to locus tag for clarity
set2genes <-
  string.enrichment |>
  select(
    database = category,
    geneset = `term name`,
    geneset.name = description,
    cluster,
    genes
  ) |>
  unique() |>
  separate_rows(genes, sep = '\\|') |>
  left_join(mcl |> select(- cluster), c('genes' = 'string.name'))

# save for later use
enrich |>
  write_tsv('analysis/L_string-enrichment.tsv')

set2genes |>
  write_tsv('analysis/L_string-enrichment-genes.tsv')

################################################################################
################################################################################
# Identify potential "correlation" of enriched genesets to expression levels
# with heatmaps

# per gene z-scale
z.vst <-
  vst.mat |>
  apply(1, scale) |>
  t() |>
  magrittr::set_colnames(colnames(vst.mat)) 


# select database of interest
sel.db <-
  set2genes |>
  filter(database == 'STRING Clusters')

# Reduce redundancy
# Remove geneset if there are larger ones of same name
# (stringApp internally uses the Jaccard index, but not relevant here)
mask <-
  sel.db |>
  select(geneset, geneset.name, Geneid) |>
  unique() |>
  count(geneset, geneset.name) %>%
  left_join(., ., 'geneset.name') |>
  filter(n.x < n.y) |>
  select(geneset = geneset.x) |>
  unique()

sel <-
  sel.db |>
  anti_join(mask, 'geneset') |>
  transmute(
    set = sprintf(
      '%s (%s)',
      str_replace(geneset.name, '(?<=.{40})....*$', '...'),
      geneset
    ),
    Geneid
  ) |>
  unique()

################################################################################
# Build heatmap annotation colors

# color CO2 condition
col.annot <-
  meta |>
  with(data.frame('CO2' = as.character(condition), row.names = lib))

# color cluster
row.annot <-
  sel |>
  select(Geneid) |>
  unique() |>
  left_join(mcl, 'Geneid') |>
  # cluster as string with leading space
  mutate(cluster = sprintf('%2s', cluster)) |>
  select(cluster, Geneid) |>
  arrange(cluster) |>
  with(data.frame('Cluster' = cluster, row.names = Geneid))

annot.cl <- list(
  'CO2' =  meta |>
    select(CO2 = condition) |>
    unique() |>
    arrange(CO2) |>
    mutate(cl = cbPalette[-5][1:n()]) |>
    with(set_names(cl, as.character(CO2)))
)

# cluster cols in ascending CO2 order
lib.clust <-
  vst.mat |>
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

################################################################################

dir.create('analysis/L_geneset-heatmaps')

helper <- function(i) {
  # i <- 'Translation, and Ribosome biogenesis (CL:2400)'
  path <- sprintf('analysis/L_geneset-heatmaps/%s.jpeg',
                  str_replace_all(i, '[^a-zA-Z0-9]', '_'))
  mask <-
    sel |>
    filter(set == i) |>
    pull(Geneid)
  
  pheatmap::pheatmap(
    vst.mat[mask, ],
    scale = 'row',
    cluster_cols = lib.clust,
    show_rownames = FALSE,
    show_colnames = FALSE,
    annotation_col = col.annot,
    annotation_row = row.annot,
    annotation_colors = annot.cl,
    color = colorRampPalette(rev(
      RColorBrewer::brewer.pal(n = 7, name = "RdBu")))(59),
    filename = path
  )
}

sel |>
  pull(set) |>
  unique() |>
  future_map(helper)

################################################################################

