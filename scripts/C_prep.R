#!/usr/bin/env Rscript

# Prepare raw input into nicer tables

library(tidyverse)

################################################################################
# Process raw counts of CO2 data


raw.counts <-
  'raw-data/co2-adaptation-counts.tsv' |>
  read_tsv() |>
  select(-c(Chr, Start, End, Strand, Length))

meta <-
  tibble(lib = colnames(raw.counts)[-1]) |>
  separate(lib, c('CO2', 'sample', 'dataset'),
           sep = '_', remove = FALSE) |>
  select(- dataset) |>
  mutate(
    CO2 = CO2 |>
      str_remove('^CO2-') |>
      str_remove('percent$') |>
      str_replace('-', '.') |>
      as.numeric(),
    lane = str_remove(sample, '-.*$'),
    sample = str_remove(sample, '^Lane.-')
  )

write_tsv(meta, 'data/C_meta.tsv')
write_tsv(raw.counts, 'data/C_raw-counts.tsv')

################################################################################
################################################################################
# Load gene annotation

gff <- 'raw-data/PCC7002-genome.gff.gz' |>
  rtracklayer::import.gff3()


gff |>
  plyranges::filter(type == 'gene') |>
  as_tibble() |>
  select(Geneid = ID, type = gene_biotype,
         locus_tag, old_locus_tag,
         name = Name) |>
  left_join(
    gff |>
      select(ID, Parent, product) |>
      as_tibble() |>
      select(Parent, product) |>
      drop_na(product) |>
      unnest(Parent) |>
      unique(),
    c('Geneid' = 'Parent')
  ) -> annot.gff

# Add Rfam hits info
raw.counts |>
  select(Geneid) |>
  filter(!str_detect(Geneid, '^gene')) |>
  mutate(
    type = 'ncRNA',
    name = str_remove(Geneid, '^Candidate_')
  ) |>
  bind_rows(annot.gff) -> annot

write_tsv(annot, 'data/C_annotation.tsv')
################################################################################
sessionInfo()
