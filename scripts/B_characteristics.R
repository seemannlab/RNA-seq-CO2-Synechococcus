# Extract wet-lab measured characteristics data from Excel sheets
# and make prettier plots

library(tidyverse)
library(ggpubr)

cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442",
               "#0072B2", "#D55E00", "#CC79A7")

set.seed(123)

################################################################################
# Extract from excel


readxl::read_excel(
  'raw-data/RNA_JAN23.xlsx',
  # only read 2 x 4 data rows plus header without summary
  n_max = 16 + 1
) -> dat

dat |>
  rename(
    "Ph" = `Ph (at collection time)`
  ) |>
  select(- contains('TIC')) |>
  rename(
    'phycocyanin (PC) / biomass' = `PC/biomass`,
    'chlorophyll a (Chl a) / biomass' = `Chl a/biomass`
  ) |>
  mutate(
    Photons = ifelse(
      `CO2%` %in% c(1, 15),
      150,
      250
    )
  ) -> data

write_tsv(data, 'data/B_characteristics.tsv')

################################################################################
# Pretty plot of growth curve


data %>%
  pull(`CO2%`) %>%
  unique %>%
  sort() %>%
  as.character() %>%
  tibble(a = ., b = lead(.)) %>%
  drop_na() %>%
  add_row(a = '0.04', b = '30') |>
  rowwise() %>%
  do(i = unlist(.)) %>%
  pull(i) -> cmps

# order for nicer overview
cmps <- cmps[c(1, 3, 2, 4)]

data  %>%
  rename(CO2 = `CO2%`) %>%
  ggboxplot(
    x = 'CO2',
    y = 'Growth rate',
    add = 'jitter', add.params = list(color = 'blue')
  ) +
  xlab('CO2 aeration %') +
  # ylab('measurement metric (different in each facet)') +
  ylab('Growth rate') +
  stat_compare_means(
    comparisons = cmps,
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

ggsave('data/B_growth.jpeg', width = 8, height = 7, dpi = 400)

################################################################################

sessionInfo()