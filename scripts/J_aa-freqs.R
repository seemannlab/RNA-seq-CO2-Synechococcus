# Investigate the impact of differential expression on amino acid compositions

library(tidyverse)
library(ggpubr)
library(patchwork)

# https://riptutorial.com/r/example/28354/colorblind-friendly-palettes
cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442",
               "#0072B2", "#D55E00", "#CC79A7")

################################################################################
################################################################################
# Load data

annot <-
  'data/C_annotation.tsv' |>
  read_tsv()

deg <-
  'analysis/D_stagewise-adjusted-DEGs.tsv' |>
  read_tsv()

norm <-
  'analysis/D_normalized-counts.tsv' |>
  read_tsv()

vst <-
  'analysis/D_vst-expression.tsv' |>
  read_tsv()

meta <-
  'data/C_meta.tsv' |>
  read_tsv()

freqs <-
  'analysis/G_frequencies.tsv' |>
  read_tsv()

freqs.mat <-
  freqs |>
  select(- Geneid) |>
  as.matrix() |>
  magrittr::set_rownames(freqs$Geneid)

################################################################################
# Overall AA freq

freqs |>
  pivot_longer(- Geneid) |>
  mutate(AA = fct_reorder(name, value)) -> foo
foo |>
  ggplot(aes(AA, value, fill = AA)) +
  geom_violin() +
  # geom_boxplot(fill = 'white', width = .2) +
  xlab(NULL) +
  ylab('Gene length normalized\nAmino acid frequency') +
  theme_pubr(18) +
  theme(
    legend.position = 'hide',
    axis.text.x = element_text(angle = 60, hjust = 1)
  ) -> p1.total


################################################################################
# Z-scale frequency per amino acid
fz <-
  freqs.mat |>
  apply(2, scale) |>
  magrittr::set_rownames(rownames(freqs.mat))

fz |>
  as_tibble(rownames = 'Geneid') |>
  pivot_longer(- Geneid) |>
  mutate(AA = fct_relevel(name, levels(foo$AA))) |>
  ggplot(aes(AA, value, fill = AA)) +
  geom_boxplot() +
  xlab(NULL) +
  ylab('z-scaled gene length normalized\nAmino acid frequency') +
  theme_pubr(18) +
  theme(
    legend.position = 'hide',
    axis.text.x = element_text(angle = 60, hjust = 1)
  ) -> p1.z

p1.total | p1.z

ggsave('analysis/J_freqs-overall.jpeg', width = 10, height = 7, dpi = 400)

################################################################################
# Explore distribution of AA frequencies


fz.long <-
  fz |>
  as_tibble(rownames = 'Geneid') |>
  pivot_longer(- Geneid, names_to = 'aa')

fz.long |>
  ggecdf('value', color = 'aa')


p3 <-
  ( p3.hist / p3.qq / p3.cor ) +
  plot_layout(heights = c(1, 1, 2)) +
  plot_annotation(tag_levels = 'A')

ggsave('analysis/G_AA-dist-cor.jpeg', p3,
       width = 12, height = 12, dpi = 400)

################################################################################
# Correlate AA and expression

# Per gene z-Scaled expression
vz <-
  vst2 |>
  apply(1, scale) |>
  t() |>
  magrittr::set_colnames(colnames(vst2))

# # keep only top X most variant expressed genes among the DEGs
# v.dat <-
#   tibble(
#     gene = rownames(vst2),
#     v =  MatrixGenerics::rowVars(vst2)
#   )
# 
# qv <- quantile(v.dat$v, .9)
# nice <-
#   sprintf(
#     'Top 10%% (n = %g genes) with\nvariance ≥ %.2f',
#     v.dat |>
#       filter(v >= qv) |>
#       nrow(),
#     qv
#   )
# 
# v.dat |>
#   ggplot(aes( v)) +
#   stat_ecdf() +
#   geom_hline(yintercept = .9, linetype = 'dotted') +
#   geom_vline(xintercept = qv, color = 'red') +
#   scale_y_continuous(breaks = seq(0, 1, .1)) +
#   annotate(
#     'text', label = nice,
#     x = 1, y = .84, color = 'red',
#     size = 8,
#     hjust = 0
#   ) +
#   xlab('Expression variance (vst)') +
#   ylab('Empirical cumulative density') +
#   theme_pubr(18)
# 
# ggsave('analysis/I_top-var.jpeg',
#        width = 8, height = 8, dpi = 400)
# 
# 
# top.var <-
#   v.dat |>
#   filter(v >= qv) |>
#   pull(gene)

# nicer library names
look <-
  meta |>
  mutate(nice = sprintf('%s%% (%s)', CO2, sample)) |>
  with(set_names(nice, lib))

deg |>
  filter(is.de, abs(log2FoldChange) >= 1) |>
  pull(Geneid) |>
  unique() |>
  intersect(rownames(f2z)) -> mask

crossing(
  AA = colnames(f2z),
  lib = colnames(vst2)
) |>
  group_by(AA, lib) |>
  # do(r2 = cor.test(f2z[top.var, .$AA], vz[top.var, .$lib])$estimate) |>
  # do(r2 = cor.test(f2z[, .$AA], vz[, .$lib])$estimate) |>
  do(r2 = cor.test(f2z[mask, .$AA], vst2[mask, .$lib])$estimate) |>
  # do(r2 = cor.test(sqrt(freqs2)[mask, .$AA], vz[mask, .$lib])$estimate) |>
  ungroup() |>
  unnest(r2) |>
  mutate(lib = look[lib]) |>
  spread(lib, r2) -> foo

foo |>
  select(- AA) |>
  as.matrix() |>
  magrittr::set_rownames(foo$AA) |>
  pheatmap::pheatmap(
    scale = 'none',
    show_colnames = TRUE,
    show_rownames = TRUE,
    display_numbers = TRUE,
    fontsize = 18,
    number_color = 'black',
    color = colorRampPalette(rev(
      RColorBrewer::brewer.pal(n = 5, name = "RdBu")))(59),
    # filename = 'analysis/I_AA-expr-cor.jpeg', width = 14, height = 14
  )
dev.off()
  


################################################################################
# Signal in frequencies

pca <- prcomp(t(f2z))
percentVar <- ( pca$sdev^2/sum(pca$sdev^2) ) %>%
  set_names(paste0('PC', 1:length(.)))

x <- 'PC1'
y <- 'PC2'
pca$x |>
  as_tibble(rownames = 'amino') |>
  left_join(aaMap, c('amino' = 'name')) |>
  ggplot(aes(!! sym(x), !! sym(y), color = scProp)) +
  # ggplot(aes(!! sym(x), !! sym(y))) +
  geom_point(size = 5) +
  ggsci::scale_color_jama(name = 'side chain property at pH 7') +
  xlab(sprintf('%s: %.1f%% variance', x, percentVar[[x]] * 100)) +
  ylab(sprintf('%s: %.1f%% variance', y, percentVar[[y]] * 100)) +
  ggrepel::geom_label_repel(aes(label = amino), show.legend = FALSE, size = 5) +
  coord_cartesian() +
  theme_pubr(18) +
  ggtitle('PCA of amino acid frequencies alone')

ggsave('analysis/I_AA-freq-pca.jpeg', width = 10, height = 8, dpi = 400)

################################################################################
################################################################################
