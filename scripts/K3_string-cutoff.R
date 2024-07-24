# Determine a cutoff for association scores

library(tidyverse)
library(ggpubr)

string <-
  'https://stringdb-downloads.org/download/protein.links.detailed.v12.0/32049.protein.links.detailed.v12.0.txt.gz' |>
  read_delim(delim = ' ') |>
  mutate_if(is.numeric, ~ .x / 1e3)

string |>
  pull(combined_score) |>
  quantile(.8) -> xs

string |>
  ggplot(aes(combined_score)) +
  scale_x_continuous(breaks = seq(0, 1, .1)) +
  scale_y_continuous(breaks = seq(0, 1, .1)) +
  geom_hline(yintercept = .8, color = 'red') +
  geom_vline(xintercept = xs, color = 'red') +
  annotate('text', x = .75, y = .1, label = xs, color = 'red', size = 5) +
  stat_ecdf() +
  xlab('STRING association score') +
  ylab('Empirical cumulative density') +
  theme_pubr(18)

ggsave('analysis/K3_string-cutoff.jpeg', width = 7, height = 6, dpi = 400)