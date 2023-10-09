#!/usr/bin/env Rscript

# Compute benchmark of the various networks
library(tidyverse)
library(ggpubr)
library(patchwork)


# https://riptutorial.com/r/example/28354/colorblind-friendly-palettes
cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442",
               "#0072B2", "#D55E00", "#CC79A7")

################################################################################

enrichment.data <- read_tsv('analysis/H_enrichment-data.tsv.gz')

################################################################################
# Plot enrichment of recall curves

enrichment.data |>
  drop_na() |>
  mutate_at('nice.name', str_replace, ', ', '\n') |>
  ggplot(aes(threshold, ratio, col = channel)) +
  geom_line(size = 1.2, alpha = 0.7) +
  scale_color_manual(values = RColorBrewer::brewer.pal(10, 'Paired')) +
  scale_x_continuous(breaks = seq(0, 1, .2)) +
  xlab('Cut-off') +
  ylab('Observed over expected recall rate') +
  # facet_grid(reference ~ network, scales = 'free') +
  facet_grid(net ~ nice.name, scales = 'free') +
  theme_bw(18) +
  theme(legend.position = 'bottom')

ggsave('analysis/H_enrich-curve.jpeg', width = 18, height = 8)

################################################################################
# Approximate ROC curves

enrichment.data |>
  mutate(
    tp = expected,
    fp = observed - expected,
    # n = pairs.total - expected.total,
    tn = pairs.total - expected.total - fp,
    fn = expected.total - tp,
    tpr = tp / (tp + fn),
    fpr = fp / (fp + tn)
  ) -> approx.roc

approx.roc |>
  select(net, nice.name, channel, fpr, tpr) |>
  mutate_at('nice.name', str_replace, ', ', '\n') |>
  ggplot(aes(fpr, tpr, color = channel)) +
  geom_line(size = 1.2, alpha = 0.7) +
  scale_x_continuous(breaks = seq(0, 1, .2)) +
  scale_y_continuous(breaks = seq(0, 1, .2)) +
  scale_color_manual(values = RColorBrewer::brewer.pal(10, 'Paired')) +
  geom_abline(slope = 1, color = 'black', linetype = 'dashed') +
  xlab('False positive rate') +
  ylab('True positive rate') +
  facet_grid(net ~ nice.name) +
  theme_bw(18) +
  theme(legend.position = 'bottom')

ggsave('analysis/H_roc-curve.jpeg', width = 18, height = 8)

################################################################################
# Rank by approximate Area under the ROC curve

approx.roc.area <-
  approx.roc |>
  select(net, nice.name, channel, fpr, tpr) |>
  arrange(net, nice.name, channel, fpr) |>
  unique() |>
  group_by(net, nice.name, channel) |>
  mutate(
    next.fpr = lead(fpr),
    next.tpr = lead(tpr),
    # difference in x
    xi = next.fpr - fpr,
    # min tpr
    tpr.min = pmin(tpr, next.tpr),
    # the height of the "triangle part"
    tpr.tri = abs(tpr - next.tpr),
    # area is about: rectangle + triangle
    area.part = xi * (tpr.min + tpr.tri / 2)
  ) |>
  summarize(roc.area = sum(area.part, na.rm = TRUE)) |>
  ungroup() |>
  arrange(roc.area)
  
write_tsv(approx.roc.area, 'analysis/H_roc-areas.tsv')

approx.roc.area |>
  spread(channel, roc.area) |>
  arrange(desc(combined_score)) |>
  write_tsv('analysis/H_roc-areas-wide.tsv')

################################################################################
# Nice visualization for the ROC areas


approx.roc.area |>
  spread(net, roc.area) |>
  mutate_at('nice.name', str_replace, ', ', '\n') |>
  ggplot(aes(`Pearson correlation`, `WGCNA TOM`, color = channel)) +
  geom_abline(slope = 1, color = 'black', linetype = 'dashed') +
  geom_point(size = 5) +
  scale_x_continuous(breaks = seq(.4, .7, .1)) +
  scale_y_continuous(breaks = seq(.4, .7, .1)) +
  coord_fixed() +
  facet_wrap(~ nice.name) +
  scale_color_manual(values = RColorBrewer::brewer.pal(10, 'Paired')) +
  # scale_color_manual(values = cbPalette) +
  theme_pubr(18) +
  ggtitle('Comparison of areas under the ROC curves')

ggsave('analysis/H_aroc-scatter.jpeg', width = 9, height = 10)
