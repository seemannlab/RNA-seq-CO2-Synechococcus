#!/usr/bin/env Rscript

# Extract wet-lab measured characteristics data from Excel sheets
# and make prettier plots

library(tidyverse)
library(ggpubr)

################################################################################
# Extract from excel

readxl::read_excel(
  'raw-data/RNA_July22.xlsx',
  # only read 2 x 4 data rows plus header without summary
  n_max = 8 + 1
) -> data1

readxl::read_excel(
  'raw-data/RNA_JAN23.xlsx',
  # only read 2 x 4 data rows plus header without summary
  n_max = 16 + 1
) -> data2

bind_rows(
  data1 %>%
    mutate_at('Sample', str_replace, 'S', '7002-') %>%
    rename(
      "Ph" = `Ph (at collection time)`
    ),
  data2 %>%
    rename(
      "OD (Abs 730)" = `OD (Abs 730) when harvest`,
      "Ph" = `Ph (at collection time)`
    )
) %>%
  select(- contains('TIC')) %>%
  rename(
    'phycocyanin (PC) / biomass' = `PC/biomass`,
    'chlorophyll a (Chl a) / biomass' = `Chl a/biomass`
  ) -> data

write_tsv(data, 'data/B_characteristics.tsv')

################################################################################
# Pretty plot

data %>%
  pull(`CO2%`) %>%
  unique %>%
  sort() %>%
  as.character() %>%
  tibble(a = ., b = lead(.)) %>%
  drop_na() %>%
  rowwise() %>%
  do(i = unlist(.)) %>%
  pull(i) -> cmps


data  %>%
  select(Sample, `CO2%`,
         `chlorophyll a (Chl a) / biomass`, `phycocyanin (PC) / biomass`,
         `Growth rate`) %>%
  # select(- `Doubling time (h)`) %>%
  gather('metric', 'value', - c( Sample, `CO2%`))  %>%
  rename(CO2 = `CO2%`) %>%
  ggboxplot(
    x = 'CO2',
    y = 'value'
    # add = 'jitter'
  ) +
  xlab('CO2 aeration %') +
  # ylab('measurement metric (different in each facet)') +
  ylab(NULL) +
  stat_compare_means(
    comparisons = cmps,
    label = 'p.signif',
    method = 't.test',
    size = 5
    
  ) + 
  facet_wrap(~ metric, scales = 'free') +
  theme_pubr(18) -> p

annotate_figure(
  p,
  top = text_grob('Experimental characteristics',
                  face = "bold", size = 18),
  bottom = text_grob(
    paste(
      'T-test P-value:',
      'ns  not significant;',
      '*  < 0.05;',
      '**  < 0.01;',
      '***  < 0.001;',
      '****  < 0.0001',
      sep = ' '
    ),
    hjust = 1, x = .9,
    face = "italic", size = 14,
  )
)

ggsave('data/B_characteristics.jpeg', width = 12, height = 8, dpi = 400)

################################################################################
# Repeat but with focus on the growth curve

# order for nicer overview
cmps2 <- cmps[c(5:3, 1:2)]

data  %>%
  rename(CO2 = `CO2%`) %>%
  ggboxplot(
    x = 'CO2',
    y = 'Growth rate'
    # add = 'jitter'
  ) +
  xlab('CO2 aeration %') +
  # ylab('measurement metric (different in each facet)') +
  ylab('Growth rate') +
  stat_compare_means(
    comparisons = cmps2,
    label = 'p.signif',
    method = 't.test',
    size = 5
    
  ) + 
  theme_pubr(18) -> p

annotate_figure(
  p,
  bottom = text_grob(
    paste(
      'T-test P-value:',
      'ns  not significant;',
      '*  < 0.05;',
      '**  < 0.01;',
      '***  < 0.001;',
      '****  < 0.0001',
      sep = ' '
    ),
    hjust = 1, x = .9,
    face = "italic", size = 14,
  )
)

ggsave('data/B_growth.jpeg', width = 8, height = 5, dpi = 400)

################################################################################