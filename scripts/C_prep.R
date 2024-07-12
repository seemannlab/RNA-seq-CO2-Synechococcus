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
  colnames(raw.counts)[-1] |>
  tibble(lib = _) |>
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
         name = Name,
         seqnames, start, end, strand) |>
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

'raw-data/co2-adaptation-counts.tsv' |>
  read_tsv() |>
  select(
    Geneid, 
    seqnames = Chr,
    start = Start, 
    end = End,
    strand = Strand
  ) |>
  filter(!str_detect(Geneid, '^gene')) |>
  mutate(
    type = ifelse(str_detect(Geneid, 'motif'), 'ncRNA-CRS', 'ncRNA-Rfam-hit'),
    name = str_remove(Geneid, '^Candidate_')
  ) |>
  bind_rows(annot.gff) -> annot

# improve type names
annot <-
  annot |>
  mutate(type = case_when(
    type == 'protein_coding' ~ 'protein-coding',
    type == 'tRNA' ~ 'tRNA',
    type == 'rRNA' ~ 'rRNA',
    type == 'ncRNA-CRS' ~ 'CRS',
    # Some manual coding
    Geneid == 'Candidate_RF01343_CRISPR-DR33' ~ 'CRISPR-DR33',
    Geneid == 'Candidate_RF01701_Cyano-1' ~ 'Cyano-1',
    Geneid == 'Candidate_RF02514_5_ureB_sRNA' ~ '5\' ureB sRNA',
    Geneid == 'Candidate_RF02701_PsrR1' ~ 'PsrR1',
    Geneid == 'gene-SYNPCC7002_RS16250' ~ 'tmRNA',
    Geneid == 'gene-SYNPCC7002_RS16255' ~ 'RNase P',
    Geneid == 'gene-SYNPCC7002_RS16460' ~ 'SRP sRNA',
    Geneid == 'gene-SYNPCC7002_RS16495' ~ '6S RNA'
  ))
# annot |> filter(type != 'protein-coding', type != 'tRNA', type != 'rRNA') |> View()

write_tsv(annot, 'data/C_annotation.tsv')
################################################################################
sessionInfo()
