#!/usr/bin/env Rscript

# Prepare raw input into nicer tables

library(tidyverse)

################################################################################
# Process raw counts of CO2 data

'raw-data/schlange-data1.csv' %>%
  read_csv() %>%
  select(sample, batch) %>%
  unique -> batches

'raw-data/1-15pct-co2-counts.tsv' %>%
  read_tsv() -> x1
'raw-data/co2-adaptation-counts.tsv' %>%
  read_tsv() -> x2

list(x1, x2) %>%
  map(select, -c(Chr, Start, End, Strand, Length)) %>%
  reduce(full_join, by = 'Geneid') -> raw.counts

tibble(lib = colnames(raw.counts)[-1]) %>%
  separate(lib, c('condition', 'sample', 'dataset'),
           sep = '_', remove = FALSE) %>%
  select(- dataset) %>%
  mutate_at('condition', str_remove, '^CO2-') %>%
  left_join(batches, 'sample') %>%
  mutate(
    batch = ifelse(
      is.na(batch),
      paste0('batch3-', str_extract(sample, 'Lane.')),
      batch
    )
  ) %>%
  mutate_at('sample', str_remove, '^Lane.-') %>%
  mutate_at('condition', str_remove, 'percent$') %>%
  mutate_at('condition', str_replace, '-', '.') %>%
  # as factors in numeric order
  mutate_at('condition', ~ fct_reorder(.x, as.numeric(.x))) -> meta

write_tsv(meta, 'data/C_meta.tsv')
write_tsv(raw.counts, 'data/C_raw-counts.tsv')

################################################################################
# Process raw counts of public data

public.counts <-
  'raw-data/public-datasets/public-raw-counts.tsv' |>
  read_tsv() |>
  select(-c(Chr, Start, End, Strand, Length))

meta.public <-
  tibble(lib = colnames(public.counts)[-1]) %>%
  separate(lib, c('condition', 'sample', 'dataset'),
           sep = '_', remove = FALSE) %>%
  mutate_at('condition', str_remove, '^replicate.-')

write_tsv(meta.public, 'data/C_public-meta.tsv')
write_tsv(public.counts, 'data/C_public-raw-counts.tsv')

################################################################################
# Load gene annotation

gff <- 'raw-data/PCC7002-genome.gff.gz' %>%
  rtracklayer::import.gff3()


gff %>%
  plyranges::filter(type == 'gene') %>%
  as_tibble() %>%
  select(Geneid = ID, type = gene_biotype, locus_tag, old_locus_tag, name = Name) %>%
  left_join(
    gff %>%
      select(ID, Parent, product) %>%
      as_tibble() %>%
      select(Parent, product) %>%
      drop_na(product) %>%
      unnest(Parent) %>%
      unique,
    c('Geneid' = 'Parent')
  ) -> annot.gff

# Add Rfam hits info
raw.counts %>%
  select(Geneid) %>%
  filter(!str_detect(Geneid, '^gene')) %>%
  mutate(
    type = 'ncRNA',
    name = str_remove(Geneid, '^Candidate_')
  ) %>%
  bind_rows(annot.gff) -> annot

write_tsv(annot, 'data/C_annotation.tsv')
################################################################################