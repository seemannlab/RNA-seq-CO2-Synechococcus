## NOTE: Cytoscape with the network session loaded needs to be running for this script

# Script is nearly identical to F_cluster-overview.R
# But run in combination with the K_30 StringApp results



library(tidyverse)

library(RCy3)

library(conflicted)

conflicts_prefer(purrr::map)
conflicts_prefer(dplyr::select)
conflicts_prefer(dplyr::filter)
conflicts_prefer(dplyr::count)

################################################################################

annot <-
  'data/C_annotation.tsv' |>
  read_tsv()

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
  write_tsv('analysis/K2_mcl-clustering.tsv')

################################################################################
# use command file to export all the clusters enrichment from cytoscape

x <- 1:10

# Query cytoscape for enrichment per cluster
c(
  'string retrieve group-wise enrichment maxGroupNumber=100 minGroupSize=5',
  sprintf(
    'table export outputFile="%s/analysis/K2_cluster-enrichment/cluster_%g.csv" table="STRING Enrichment: __mclCluster %g"',
    getwd(), x, x
  )
) |>
  write_lines('analysis/K2_cluster-enrichment-export-commands.txt')

dir.create('analysis/K2_cluster-enrichment/')
tryCatch({
  RCy3::commandRunFile(file.path(
    getwd(),
    'analysis/K2_cluster-enrichment-export-commands.txt'
  ))
},
# ignore error on "empyt tables'
error=function(error_message) {}
)

# load the exported data
string.enrichment <-
  'analysis/K2_cluster-enrichment/*csv' |>
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
  write_tsv('analysis/K2_string-enrichment.tsv')

set2genes |>
  write_tsv('analysis/K2_string-enrichment-genes.tsv')

################################################################################
################################################################################

sessionInfo()