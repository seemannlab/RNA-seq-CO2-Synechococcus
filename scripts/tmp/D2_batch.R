#!/usr/bin/env Rscript

# Inspect to what extend the batch is a concern in the PCA plots with/without air
library(tidyverse)
library(ggpubr)

library(corrplot)
library(broom)

################################################################################
# Load data

meta <-
  'data/C_meta.tsv' |>
  read_tsv()

experiment <-
  'data/B_characteristics.tsv' |>
  read_tsv()


################################################################################
# Expression matrix

foo <-
  'analysis/E_vst.tsv' |>
  read_tsv()

dba <-
  foo |>
  select(- Geneid) |>
  as.matrix() |>
  magrittr::set_rownames(foo$Geneid)

################################################################################
# Data without air

meta2 <-
  meta |>
  filter(condition != 0.04)
dba2 <- dba[, meta2$lib]

################################################################################
# Run PCA

helper <- function(dba, meta, ntop = 500) {
  # Select most variant genes
  rv <- MatrixGenerics::rowVars(dba)
  sel <- order(rv, decreasing = TRUE)[seq_len(min(ntop, length(rv)))]
  # Compute PCA
  pca <- prcomp(t(dba[sel, ]))

  # Table of PCA and experimental data
  pca$x |>
    as_tibble(rownames = 'lib') |>
    select(lib, PC1, PC2, PC3, PC4) |>
    left_join(
      meta |>
        select(lib, batch, sample),
      'lib'
    ) |>
    left_join(
      transmute(
        experiment,
        sample = Sample,
        CO2 = `CO2%`,
        `Growth rate`,
        # 'Chlorophyll A' = `chlorophyll a (Chl a) / biomass`,
        # Phycocyanin = `phycocyanin (PC) / biomass`
      ),
      'sample'
    ) |>
    select(- lib, - sample) |>
    mutate(
      CO2 = CO2 |>
        as_factor(),
      batch = batch |>
        str_remove('-Lane.$') |>
        as_factor()
    )
}

dat1 <- helper(dba, meta)
dat2 <- helper(dba2, meta2)

################################################################################
# Pairwise correlations

c1 <-
  dat1 |>
  mutate(CO2 = CO2 |>
           as.character() |>
           as.numeric()) |>
  select_if(is.numeric) |>
  as.matrix() |>
  cor()

c2 <-
  dat2 |>
  mutate(CO2 = CO2 |>
           as.character() |>
           as.numeric()) |>
  select_if(is.numeric) |>
  as.matrix() |>
  cor()

# Combine in one matrix
c1[lower.tri(c1)] <- c2[lower.tri(c2)]

jpeg('analysis/D2_correlations.jpeg', res = 400,
     width = 20, height = 20, units = 'cm')
corrplot::corrplot(
  c1,
  method = 'square',
  order = 'original',
  addCoef.col = 'black',
  mar = c(0, 0, 3, 0),
  col = colorRampPalette(c('#FF6666', '#DDDDDD', '#6666FF'))(50),
  title = 'Correlations when including air (upper) vs excluding (lower triangle)'
)
dev.off()

################################################################################
# Change in residuals when excl variables

## Code below dynamically creates something like code below but for all PCs and
## both dat1 and dat2
# list(
#   full   = glm(PC1 ~ batch + CO2 + `Growth rate` , data = dat),
#   batch  = glm(PC1 ~         CO2 + `Growth rate` , data = dat),
#   co2    = glm(PC1 ~ batch +       `Growth rate` , data = dat),
#   growth = glm(PC1 ~ batch + CO2 +               , data = dat)
# ) |>
#   map(~ sum(.x$residuals ** 2)) |>
#   as_tibble() |>
#   mutate(PC = 'PC1')

var_combinations <- function(dat, dataset, pcx) {
  # pcx <- 'PC1'
  vars <- list('batch', 'CO2', '`Growth rate`')
  
  # All combinations of variables
  combs <-
    vars |>
    seq_along() |>
    # extract combinations from the columns
    map(~ combn(vars, .x)) |>
    map(apply, 2, c) |>
    reduce(c) |>
    map(unlist)
  
  # names from variables
  ns <-
    combs |>
    map(str_remove_all, '`') |>
    map(str_remove_all, ' .*$') |>
    map(str_c, collapse = '+')

  # Dynamically build formulas
  fmls <- 
    combs |>
    map(~ paste(pcx, '~', str_c(.x, collapse = ' + '))) |>
    set_names(ns)
  
  # Run GLM on dataset
  fmls |>
    map(glm, data = dat) |>
    # Collect sum of square residuals
    # map(~ sum(.x$residuals ** 2)) |>
    map('deviance') |>
    as_tibble() |>
    # Remember PC and dataset
    mutate(PC = pcx, dataset = dataset)
}

tribble(
  ~ dataset,     ~ dat,
  'Full dataset',  dat1,
  'Excluding air', dat2
) |>
  group_by_all() |>
  tidyr::expand(pcx = paste0('PC', 1:4)) |>
  ungroup() |>
  pmap(var_combinations) |>
  bind_rows() |>
  pivot_longer(- c(PC, dataset)) |>
  mutate_at('name', str_remove_all, '`') |>
  mutate(
    value = log10(value),
    name = fct_reorder(name, value, .fun = min),
    nice = sprintf('%.3f', value),
    nice.color = ifelse(value < 2, 'black', 'white')
  ) |>
  # ggplot(aes(log(value))) + geom_histogram() + facet_wrap(~ dataset)
  ggplot(aes(PC, name, fill = value, label = nice)) +
  geom_tile() +
  geom_text(aes(color = I(nice.color)), size = 5) +
  scale_fill_viridis_c(name = 'log10 sum of squared residuals', direction = -1) +
  labs(x = NULL, y = NULL) +
  ggtitle(
    'Generalized Linear Model regression',
    'Comparison for all combinations of variables'
  ) +
  facet_wrap(~ fct_rev(dataset)) +
  theme_pubclean(18)

ggsave('analysis/D2_glm-combinations.jpeg', width = 10, height = 10)

################################################################################
