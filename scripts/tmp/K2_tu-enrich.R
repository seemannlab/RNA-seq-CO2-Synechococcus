#!/usr/bin/env Rscript

# Check if neighboring genes or genes in transcriptional units are enriched
# for association scores in string

library(tidyverse)
library(ggpubr)

library(plyranges)
library(conflicted)

conflicts_prefer(dplyr::rename)
conflicts_prefer(plyranges::filter)

################################################################################
# Load data

annot <-
  'data/C_annotation.tsv' |>
  read_tsv()

string <-
  'data/H_string/associations.txt.gz' |>
  read_delim(delim = ' ') |>
  mutate_if(is.numeric, ~ .x / 1e3)

gff <- 'raw-data/PCC7002-genome.gff.gz' |>
  rtracklayer::import.gff3() |>
  plyranges::filter(type == 'gene')


################################################################################
# Parse TUs

tus <-
  'cyanocyc-transcription-units-of-Synechococcus-sp.-PCC-7002.txt' |>
  read_tsv() |>
  rename(
    tu = `Transcription-Units`,
    locus = `Genes of transcription unit`
  ) |>
  separate_rows(locus, sep = ' // ') |>
  inner_join(annot, c('locus' = 'locus_tag')) |>
  select(tu, Geneid)

################################################################################
# Find neighboring or overlapping genes

genes <-
  gff |>
  filter(gene_biotype == 'protein_coding') |>
  select(Geneid = ID) |>
  mutate(gs = strand)

tu.ranges <-
  genes |>
  as_tibble() |>
  right_join(tus, 'Geneid') |>
  group_by(seqnames, tu) |>
  summarize(
    start = min(start),
    end = max(end),
    strand = unique(strand)
  ) |>
  drop_na() |>
  as_granges()

tu.gaps <-
  plyranges::setdiff_ranges_directed(
    tu.ranges,
    GenomicRanges::reduce(genes)
  ) |>
  width()



# What are the gap sizes?
gene.gaps <-
  list(
    'strand specific' = genes,
    'strand unspecific' = mutate(genes, strand = '*')
  ) |>
  map(GenomicRanges::reduce) |>
  map(gaps) |>
  map(width) %>%
  map2(names(.), ~ tibble(mode = .y, len = .x)) |>
  bind_rows() |>
  bind_rows(tibble(mode = 'on TU', len = tu.gaps))
gene.gaps |>
  ggplot(aes(len, fill = fct_rev(mode))) +
  geom_density(alpha = .7) +
  scale_x_log10() +
  xlab('Intergenic region length') +
  annotation_logticks(sides = 'b') +
  geom_vline(xintercept = 500, color = 'red') +
  ggsci::scale_fill_jama(name = 'Gap strand') +
  theme_pubr(18)
  
ggsave('analysis/K2_gap-lengths.jpeg',
       width = 7, height = 5, dpi = 400)


neighbors <-
  c(
    flank(genes, 500, start = TRUE),
    flank(genes, 500, start = FALSE)
  ) |>
  join_overlap_intersect() |>
  filter(Geneid.x != Geneid.y)

overs <-
  join_overlap_intersect(genes, genes) |>
  filter(Geneid.x != Geneid.y)

################################################################################
# Pairwise lists with nice names

pairs <-
  list(
    'CyanoCyc TUs' =  tus %>%
      left_join(., ., 'tu', relationship = "many-to-many") |>
      filter(Geneid.x != Geneid.y) |>
      select(Geneid.x, Geneid.y),
    'Neighboring genes, same strand' = neighbors |>
      filter(gs.x == gs.y) |>
      as_tibble() |>
      select(Geneid.x, Geneid.y),
    'Neighboring genes' = neighbors |>
      as_tibble() |>
      select(Geneid.x, Geneid.y),
    'Overlapping genes, same strand' = overs |>
      filter(gs.x == gs.y) |>
      as_tibble() |>
      select(Geneid.x, Geneid.y),
    'Overlapping genes' = overs |>
      as_tibble() |>
      select(Geneid.x, Geneid.y)
  ) |>
  # make X < Y to not include symmetric cases (not counting twice)
  map(~ bind_rows(.x, .x)) |>
  map(filter, Geneid.x < Geneid.y) |>
  # Convert to string id
  map(left_join, annot, c('Geneid.x' = 'Geneid')) |>
  map(left_join, annot, c('Geneid.y' = 'Geneid')) |>
  map(select, protein1 = old_locus_tag.x, protein2 = old_locus_tag.y) |>
  map(drop_na) |>
  map(mutate_all, ~ paste0('32049.', .x))

# pairs |>
#   map(nrow) |>
#   unlist()
# CyanoCyc TUs Neighboring genes, same strand              Neighboring genes 
# 1776                           6974                          10610 
# Overlapping genes, same strand              Overlapping genes 
# 350                            582 
  
################################################################################
# Check for potential enrichment (related to H_benchmark.R)

pairs.bench <-
  bind_cols(
    # no. pairs
    pairs %>%
      map2(names(.), ~ tibble(name = .y, pairs = nrow(.x))) |>
      bind_rows(),
    # pairs from channel
    pairs |>
      map(inner_join, string, c('protein1', 'protein2')) |>
      map(select_if, is.numeric) |>
      map(map, sum) |>
      bind_rows()
  ) |> 
  pivot_longer(- c(name, pairs),
               names_to = 'channel', values_to = 'ref') |>
  left_join(
    string |>
      select_if(is.numeric) |>
      map(sum) %>%
      map2(names(.), ~ tibble(channel = .y, channel.total = .x)) |>
      bind_rows(),
    'channel'
  ) |>
  mutate(
    pairs.total = choose(length(genes), 2),
    expected = pairs / pairs.total,
    recall = ref / channel.total,
    ratio = recall / expected
  )

################################################################################
# nice summary plot

pairs.bench |>
  mutate(
    name = str_replace(name, ', ', '\n'),
    nice = sprintf('%.1f', ratio),
    cl = ifelse(ratio > 70, 'black', 'white')
  ) |>
  ggplot(aes(name, channel, fill = ratio,
             color = I(cl),
             label = nice)) +
  geom_tile() +
  scale_fill_viridis_c(name = 'Enrichment observed/expected recall') +
  guides(fill = guide_colorbar(barwidth =  unit(5, 'cm'))) +
  ylab('STRING associations') +
  xlab(NULL) +
  geom_text() +
  theme_pubr(18) +
  theme(axis.text.x = element_text(angle = 60, hjust = 1))

ggsave('analysis/K2_association-enrichment.jpeg',
       width = 10, height = 8, dpi = 400)
