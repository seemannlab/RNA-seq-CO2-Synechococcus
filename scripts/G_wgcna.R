#!/usr/bin/env Rscript

# Conduct (i) WGCNA and (ii) Pearson clustering both with/without public data
# Also, both for norm counts and vst

library(tidyverse)

library(ggpubr)
library(patchwork)

library(WGCNA)

# https://riptutorial.com/r/example/28354/colorblind-friendly-palettes
cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442",
               "#0072B2", "#D55E00", "#CC79A7")

# PRECISION / number of digits for adjacency matrix
PREC <- 3

################################################################################
# Load data

annot <-
  'data/C_annotation.tsv' |>
  read_tsv()

deg.res <-
  'analysis/E_dge-stagewise-analysis.tsv' |>
  read_tsv()

co2.vst <-
  'analysis/E_vst.tsv' |>
  read_tsv()
co2.norm <-
  'analysis/E_normalized-counts.tsv' |>
  read_tsv()
co2.meta <-
  'data/C_meta.tsv' |>
  read_tsv() |>
  mutate_at('condition', ~ fct_reorder(as.character(.x), as.numeric(.x)))

public.vst <-
  'analysis/F_public-vst.tsv' |>
  read_tsv()
public.norm <-
  'analysis/F_public-normalized-counts.tsv' |>
  read_tsv()
public.meta <-
  'data/C_public-meta.tsv' |>
  read_tsv()

################################################################################
# Concern on public data
# 7x more public libraries with unequal condition replicate numbers

# average per condition and dataset to balance a bit to a 2x more ratio

helper.avg <- function(x) {
  x |>
    pivot_longer(- Geneid, names_to = 'lib') |>
    left_join(public.meta, 'lib') |>
    mutate(grp = paste(dataset, condition, sep = ';')) |>
    group_by(Geneid, grp) |>
    summarize(value = mean(value)) |>
    ungroup() |>
    pivot_wider(names_from = 'grp')
}

public.norm.avg <-
  public.norm |>
  helper.avg()

public.vst.avg <-
  public.vst |>
  helper.avg()

################################################################################
# Expression matrices

helper.mat <- function(x) {
  x |>
    select(- Geneid) |>
    as.matrix() |>
    magrittr::set_rownames(x$Geneid)
}

dat.mat <- list(
  co2.norm = co2.norm |> mutate_if(is.numeric, log1p),
  co2.vst = co2.vst,
  public.norm = public.norm |> mutate_if(is.numeric, log1p),
  public.vst = public.vst,
  public.norm.avg = public.norm.avg |> mutate_if(is.numeric, log1p),
  public.vst.avg = public.vst.avg
) |>
  map(helper.mat)

################################################################################
# Matrices of DEGs, and "full" version that combines public and co2

des <-
  deg.res |>
  filter(is.de) |>
  pull(Geneid) |>
  unique()

dat.mat2 <-
  dat.mat |>
  map(~ .x[des, ]) |>
  # transpose for WGCNA
  map(t)

de.mat <- with(dat.mat2, list(
    'CO2, size-factor normalized' = co2.norm,
    'CO2, variance stabilized'    = co2.vst,
    'Incl public data, size-factor normalized' = rbind(co2.norm, public.norm),
    'Incl public data, variance stabilized'    = rbind(co2.vst, public.vst),
    'Incl public averaged data, size-factor normalized' = rbind(co2.norm, public.norm.avg),
    'Incl public averaged data, variance stabilized'    = rbind(co2.vst, public.vst.avg)
  ))

################################################################################
# Save Pearson correlations

helper.adjacency <- function(x, name, dir = 'G_pearson-cor') {
  path <- sprintf(
    'analysis/%s/%s.tsv.gz', 
    dir,
    str_replace_all(name, '[^A-Za-z0-9]', '-')
  )
                  
  # only export first PREC digits to reduce file size
  x %>%
    `*`(10 ** PREC) |>
    round() |>
    as_tibble(rownames = 'Geneid') |>
    write_tsv(path)
}

de.mat|>
  map(cor) %>%
  map2(names(.), ~ helper.adjacency(.x, .y))

################################################################################
################################################################################
# Determine Power

power.helper <- function(x.de, n) {
  sel <- pickSoftThreshold(
    x.de,
    nBreaks = 7, RsquaredCut = .75,
    powerVector = c(1:20)
  )
  
  sel$fitIndices %>%
    ggplot(aes(Power,  SFT.R.sq)) +
    geom_point() +
    geom_line() +
    ylab('Scale-free model R2') +
    geom_hline(yintercept = .75, color = 'blue') +
    scale_y_continuous(breaks = seq(0, 1, .1)) +
    geom_vline(xintercept = sel$powerEstimate, color = 'red') +
    theme_pubclean(18) -> p1
  
  
  sel$fitIndices %>%
    ggplot(aes(Power,  mean.k.)) +
    geom_point() +
    geom_line() +
    ylab('Average connectivity') +
    geom_vline(xintercept = sel$powerEstimate, color = 'red') +
    theme_pubclean(18) -> p2
  
  list(
    power =  sel$powerEstimate,
    p = (p1 + p2) + plot_annotation(
      title = n,
      # tag_level = 'A',
      theme = theme_pubr(18)
    )
  )
}

power.est <- 
  de.mat %>%
  map2(names(.), ~ power.helper(.x, .y))

power.est |>
  map('p') |>
  invoke(.f = cowplot::plot_grid, ncol = 2) |>
  ggsave(filename = 'analysis/G_wgcna-power-estimates.jpeg',
         dpi = 400, width = 16, height = 14)

################################################################################
# Semi-manual WGCNA run with saving the TOM adjacency

de.mat |>
  map(cor) |>
  map(abs) |>
  map2(power.est |> map('power'), ~ .x ** .y) |>
  map(function(x) {
    TOMsimilarity(x) |>
      magrittr::set_rownames(rownames(x)) |>
      magrittr::set_colnames(colnames(x))
  }) %>%
  map2(names(.), ~ helper.adjacency(.x, .y, dir = 'G_wgcna-tom'))


################################################################################
# Full WGCNA runs, incl clustering

nets <- map2(
  de.mat,
  power.est |> map('power'),
  ~ blockwiseModules(.x, power = .y)
)

################################################################################
# Save info on clustering

nets.mods <-
  nets |>
  map(function(net) {
    net$colors |>
      as_tibble(rownames = 'Geneid') |>
      mutate_at('value', fct_infreq) |>
      arrange(value) %>%
      left_join(
        . |>
          count(value)  |>
          mutate(m = 1:n()),
        'value'
      ) |>
      mutate(module = sprintf('WGCNA-Module:%2g', m)) |>
      select(Geneid, module)
  }) 


modules <-
  nets.mods %>%
  map2(names(.), ~ mutate(.x, net = .y)) |>
  bind_rows()

write_tsv(modules, 'analysis/G_wgca-modules.tsv')
            
################################################################################

modules |>
  select(net, module) |>
  unique() |>
  count(net, name = 'No. modules')
# net                                               `No. modules`
# 1 CO2, size-factor normalized                                  27
# 2 CO2, variance stabilized                                     32
# 3 Incl public averaged data, size-factor normalized             3
# 4 Incl public averaged data, variance stabilized                3
# 5 Incl public data, size-factor normalized                      2
# 6 Incl public data, variance stabilized                         2
