#!/usr/bin/env Rscript

# Compute benchmark of the various networks
library(tidyverse)
library(ggpubr)
library(patchwork)

library(furrr)
plan(multisession)

# https://riptutorial.com/r/example/28354/colorblind-friendly-palettes
cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442",
               "#0072B2", "#D55E00", "#CC79A7")

# PRECISION / number of digits for adjacency matrix
# (should be same as in G_wgcna.R)
PREC <- 3

################################################################################

string.dat <-
  'data/H_string/associations.txt.gz' |>
  read_delim(delim = ' ')

deg <-
  'analysis/E_dge-stagewise-analysis.tsv' |>
  read_tsv()

################################################################################
# match Geneid for DEGs

loci <-
  deg |>
  filter(is.de) |>
  select(Geneid, locus = old_locus_tag) |>
  drop_na() |>
  unique() |>
  mutate_at('locus', ~ paste0('32049.', .x))

string.genes <-
  string.dat |>
  inner_join(loci, c('protein1' = 'locus')) |>
  rename(gene1 = Geneid) |>
  inner_join(loci, c('protein2' = 'locus')) |>
  rename(gene2 = Geneid) |>
  select(- contains('protein')) |>
  select(gene1, gene2, everything()) |>
  mutate_if(is.numeric, ~ .x / 1000)

# string.dat |>
#   select_if(is.character) |>
#   unlist() |>
#   unique() |>
#   length()
# 2815
# string.genes |>
#   select_if(is.character) |>
#   unlist() |>
#   unique() |>
#   length()
# 2735

# % genes match-able by locus
# > 2735 / 2815 * 100
# [1] 97.15808

# % interactions match-able
# nrow(string.genes) /  nrow(string.dat) * 100
# 98.32413


################################################################################
# Build reference matrices

string.mat <-
  string.genes |>
  select_if(is.numeric) |>
  colnames() %>%
  set_names(.) |>
  map(~ select(string.genes, gene1, gene2, value := {{ .x }})) |>
  map(pivot_wider, names_from = gene2, values_fill = 0) |>
  map(function(x) {
    x |>
      select(- gene1) |>
      as.matrix() |>
      magrittr::set_rownames(x$gene1)
  })

################################################################################
# Score distributions

p <-
  string.mat %>%
  map2(names(.), ~ tibble(xi = .x[upper.tri(.x)], channel = .y)) |>
  bind_rows() |>
  filter(xi > 0) |>
  ggplot(aes(xi, fill = channel)) +
  geom_histogram() +
  xlab('STRING evidence channel score\n(exact zero omitted)') +
  scale_fill_manual(values = cbPalette) +
  facet_wrap(~ channel, scales = 'free_y') +
  theme_pubclean(18) +
  theme(legend.position = 'hide')

ggsave('analysis/H_string-hist.jpeg',
       p, dpi = 400,
       width = 14, height = 8)

################################################################################
# Load networks

# list of files to load
tasks <-
  crossing(
    dir = c('G_pearson-cor', 'G_wgcna-tom'),
    nice.name = c(
      'CO2, size-factor normalized',
      'CO2, variance stabilized',
      'Incl public data, size-factor normalized',
      'Incl public data, variance stabilized',
      'Incl public averaged data, size-factor normalized',
      'Incl public averaged data, variance stabilized'
    )
  ) |>
  mutate(
    path = sprintf(
      'analysis/%s/%s.tsv.gz', 
      dir,
      str_replace_all(nice.name, '[^A-Za-z0-9]', '-')
    ),
    net = fct_recode(
      dir,
      'Pearson correlation' = 'G_pearson-cor',
      'WGCNA TOM' = 'G_wgcna-tom'
    ),
    task.row = paste0('row', 1:n())
  )

# Load networks
tasks.mat <-
  tasks |>
  with(set_names(path, task.row)) |>
  future_map(function(p) {
    x <-
      p |>
      read_tsv()
    res <-
      x |>
      select(- Geneid) |>
      as.matrix() |>
      magrittr::set_rownames(x$Geneid)
    # scale back to ±1 scale, according to export PREC
    res <- res / (10 ** PREC)
    # keep scale 0...1 unsigned
    abs(res)
  })

# tasks.mat |>
#   map(c) |>
#   map(summary)
# confirmed 0...1 scale

################################################################################
# Investigate and transform edge weights (should not be > 1 for benchmark function!)

weight.dist <-
  tasks.mat %>%
  map2(names(.), ~ tibble(xi = .x[upper.tri(.x)], task.row = .y)) |>
  bind_rows() |>
  left_join(tasks, 'task.row')

p <-
  weight.dist |>
  mutate_at('nice.name', str_replace, ', ', '\n') |>
  ggplot(aes(xi, fill = nice.name)) +
  geom_histogram() +
  xlab('Network edge weight') +
  scale_color_manual(values = RColorBrewer::brewer.pal(10, 'Paired')) +
  #ggsci::scale_color_jco() +
  #scale_fill_manual(values = cbPalette) +
  facet_grid(net ~ nice.name, scales = 'free') +
  theme_pubclean(18) +
  theme(legend.position = 'hide')

ggsave('analysis/H_networks-hist.jpeg',
       p, dpi = 400,
       width = 21, height = 10)

################################################################################
# My Enrichment helper


my.bench <- function(row, channel) {
  # channel <- 'experimental'
  # row <- 'row1'

  # Select relevant matrices
  ref <- string.mat[[channel]]
  net <- tasks.mat[[row]]
  
  # Focus on shared genes
  shared <- intersect(colnames(ref), colnames(net))
  ref <- ref[shared, shared]
  net <- net[shared, shared]
  
  # Count recall etc per bin cutoff
  up.mask <- upper.tri(ref)
  seq(0, 1, length.out = 200) |>
    map(function(xi) {
      xi.mask <- (net >= xi) & up.mask
      tibble(
        threshold = xi,
        expected = sum(ref[xi.mask]),
        observed = sum(xi.mask)
      )
    }) %>%
    bind_rows() |>
    mutate(
      pairs.total = sum(up.mask),
      expected.total = sum(ref[up.mask]),
      shared = length(shared),
      prop.called = observed / pairs.total,
      recall = expected / expected.total,
      ratio = recall / prop.called,
      channel = channel,
      task.row = row
    )
}

################################################################################
# Run for all combinations

todo <- crossing(
  row = tasks |>
    pull(task.row),
  channel = string.mat |>
    names()
)

options(future.globals.maxSize=4000000000)
enrichment.data <-
  todo |>
  future_pmap(safely(my.bench)) |>
  map('result') |>
  bind_rows() |>
  left_join(tasks, 'task.row')

################################################################################
# save benchmark table out
# keep plotting separate, such that small changes on plotting script don't
# force a re-run of the benchmark

write_tsv(enrichment.data, 'analysis/H_enrichment-data.tsv.gz')
