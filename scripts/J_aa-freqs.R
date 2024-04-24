# Investigate the impact of differential expression on amino acid compositions

library(tidyverse)
library(ggpubr)
library(patchwork)

# library(forecast)

# https://riptutorial.com/r/example/28354/colorblind-friendly-palettes
cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442",
               "#0072B2", "#D55E00", "#CC79A7")

# nicer amino acid names
data(aaMap, package = 'Biobase')

################################################################################
################################################################################
# Load data

annot <-
  'data/C_annotation.tsv' |>
  read_tsv()

deg <-
  'analysis/D_stagewise-adjusted-DEGs.tsv' |>
  read_tsv()

vst <-
  'analysis/D_vst-expression.tsv' |>
  read_tsv()

meta <-
  'data/C_meta.tsv' |>
  read_tsv() |>
  mutate_at('CO2', ~ fct_reorder(as.character(.x), as.numeric(.x)))

freqs <-
  'analysis/G_frequencies.tsv' |>
  read_tsv()

################################################################################
# Convert to matrices

freqs.mat <-
  freqs |>
  select(- Geneid) |>
  as.matrix() |>
  magrittr::set_rownames(freqs$Geneid)

vst.mat <-
  vst |>
  select(- Geneid) |>
  as.matrix() |>
  magrittr::set_rownames(vst$Geneid)

################################################################################
# Overall AA freq

freqs |>
  pivot_longer(- Geneid) |>
  mutate(AA = fct_reorder(name, value)) -> foo
foo |>
  ggplot(aes(AA, value, fill = AA)) +
  # geom_violin() +
  # # geom_boxplot(fill = 'white', width = .2) +
  geom_boxplot() +
  xlab(NULL) +
  ylab('Gene length normalized\nAmino acid frequency') +
  theme_pubr(18) +
  theme(
    legend.position = 'hide',
    axis.text.x = element_text(angle = 60, hjust = 1)
  ) -> p1.total


################################################################################
# robust Z-scale frequency per amino acid
fz <-
  freqs.mat |>
  apply(2, function(x) {(x - median(x)) / mad(x)}) |>
  magrittr::set_rownames(rownames(freqs.mat))

fz |>
  as_tibble(rownames = 'Geneid') |>
  pivot_longer(- Geneid) |>
  mutate(AA = fct_relevel(name, levels(foo$AA))) |>
  ggplot(aes(AA, value, fill = AA)) +
  # geom_violin() +
  # geom_boxplot(fill = 'white', width = .2) +
  geom_boxplot() +
  xlab(NULL) +
  ylab('Per Amino acid robust z-scaled frequencies') +
  theme_pubr(18) +
  theme(
    legend.position = 'hide',
    axis.text.x = element_text(angle = 60, hjust = 1)
  ) -> p1.z

################################################################################
# explore transformation in combination with z-scale

# as vector
x <- c(fz) - min(fz) + 1
# transform
lmb <- BoxCox.lambda(x)
x2 <- BoxCox(x, lmb)
# convert back to matrix
ftrans <-
  x2 |>
  matrix(nrow = nrow(freqs.mat), ncol = ncol(freqs.mat)) |>
  magrittr::set_rownames(rownames(freqs.mat)) |>
  magrittr::set_colnames(colnames(freqs.mat))

ftrans <- (ftrans - median(ftrans)) / mad(ftrans)


ftrans |>
  as_tibble(rownames = 'Geneid') |>
  pivot_longer(- Geneid) |>
  mutate(AA = fct_relevel(name, levels(foo$AA))) |>
  ggplot(aes(AA, value, fill = AA)) +
  # geom_violin() +
  # geom_boxplot(fill = 'white', width = .2) +
  geom_boxplot() +
  # indicate 'standard deviations'
  geom_hline(yintercept = c(-2, 2), linetype = 'dashed') +
  xlab(NULL) +
  ylab('BoxCox transformed') +
  theme_pubr(18) +
  theme(
    legend.position = 'hide',
    axis.text.x = element_text(angle = 60, hjust = 1)
  ) -> p1.trans

# (p1.total | p1.trans | p1.z) +
(p1.total | p1.z | p1.trans) +
  plot_annotation(tag_levels = 'A')

ggsave('analysis/J_freqs-overall.jpeg', width = 17, height = 7, dpi = 400)

################################################################################
# Filter genese if AA does not deviate from average by 2 std, in any AA

mask.aa <-
  fz |>
  abs() |>
  apply( 1, max) %>%
  `>=`(2) |>
  which() |>
  names()

################################################################################
# Justify cutoff for lowely expressed, not expression changinge genes

deg |>
  filter(is.de) |>
  group_by(Geneid, name, product, baseMean) |>
  summarize(max.alf = max(abs(log2FoldChange))) |>
  ungroup() -> dat

# top x% with some rounding on cutoff
dat |>
  select_if(is.numeric) |>
  # map(quantile, probs = seq(0, 1, .1))
  map2(c(.1, .8), ~ quantile(.x, probs = .y)) |>
  map(~ set_names(.x, NULL)) |>
  unlist() |>
  round() |>
  as.list() -> xs
xs$baseMean <- round(xs$baseMean, -2) 
xs

dat |>
  ggplot(aes(baseMean, max.alf)) +
  geom_point(alpha = .5) +
  geom_hline(yintercept = xs$max.alf, color = 'red') +
  geom_vline(xintercept = xs$baseMean, color = 'blue') +
  scale_y_continuous(breaks = 0:8) +
  scale_x_log10(labels = scales::label_comma()) +
  annotate('text', x = xs$baseMean + 100, y = 0, size = 8,
           hjust = 0,
           label = xs$baseMean, color = 'blue') +
  annotate('text', x = 10, y = 0, size = 8,
           label = xs$max.alf, color = 'red') +
  xlab('Average expression (DESeq2 baseMean)') +
  ylab('Max. abs. logFC') +
  theme_pubr(18) -> p
  
p <- ggExtra::ggMarginal(p, fill = 'grey')
ggsave('analysis/J_expression-logFC.jpeg',plot = p,
       width = 9, height = 7, dpi = 400)

################################################################################
# Expression mask

mask.deg <-
  deg |>
  filter(
    is.de,
    abs(log2FoldChange) >= xs$max.alf,
    baseMean >= xs$baseMean
  ) |>
  pull(Geneid) |>
  unique()


my.mask <- intersect(mask.aa, mask.deg)

length(mask.aa)
length(mask.deg)
length(my.mask)

# > length(mask.aa)
# [1] 1961
# > length(mask.deg)
# [1] 985
# > length(my.mask)
# [1] 626

################################################################################
# Per gene z-Scaled expression

vz <-
  vst.mat |>
  apply(1, scale) |>
  t() |>
  magrittr::set_colnames(colnames(vst.mat))

################################################################################
# Data from Akashi 2002 on metabolic cost of AA
# https://www.ncbi.nlm.nih.gov/pmc/articles/PMC122586/

aa.cost <-
  tribble(
    ~ AA, ~ cost,
    'Ala', 11.7,
    'Cys', 24.7,
    'Asp', 12.7,
    'Glu', 15.3,
    'Phe', 52.0,
    'Gly', 11.7,
    'His', 38.3,
    'Ile', 32.3,
    'Lys', 30.3,
    'Leu', 27.3,
    'Met', 34.3,
    'Asn', 14.7,
    'Pro', 20.3,
    'Gln', 16.3,
    'Arg', 27.3,
    'Ser', 11.7,
    'Thr', 18.7,
    'Val', 23.3,
    'Trp', 74.3,
    'Tyr', 50.0,
  ) |>
  # rank by metabolic cost
  arrange(cost) |>
  mutate(
    cost.rank = rank(cost),
    # convert to long names
    AA = str_to_lower(AA),
    AA = with(aaMap, set_names(name, let.3))[AA]
  )

################################################################################
################################################################################
# Heatmap of AA and expression correlation

crossing(
  AA = colnames(fz),
  lib = colnames(vz)
) |>
  group_by(AA, lib) |>
  do(r2 = cor.test(fz[my.mask, .$AA], vz[my.mask, .$lib])$estimate) |>
  ungroup() |>
  unnest(r2) |>
  spread(lib, r2) -> dat

dat.cor <- dat |>
  select(- AA) |>
  as.matrix() |>
  magrittr::set_rownames(dat$AA)

# filter out low correlations
dat.cor |>
  abs() %>%
  apply(1, max) -> foo
dat.cor <- dat.cor[foo >= 0.2, ]

################################################################################
# Produce heatmap of correlations AA freq and Expression

dat.mat <- dat.cor
# rename to short sample name
look <-
  meta |>
  with(set_names(sample, lib))
colnames(dat.mat) <- look[colnames(dat.mat)]

# arrange by CO2 expression
lib.clust <-
  dat.mat |>
  t() |>
  dist() |>
  hclust() |>
  ape::as.phylo() |>
  ape::rotateConstr(
    meta |>
      arrange(CO2, sample) |>
      pull(sample)
  ) |>
  as.hclust()

cl <- list(
  'CO2' = meta |>
    pull(CO2) |>
    unique() |>
    sort() |>
    as.character() %>%
    set_names(cbPalette[c(1, 6, 2, 7)], .),
  'Energy' = c('white', 'black')
)
with(
  meta,
  data.frame('CO2' = as.character(CO2), row.names = sample)
) -> col.df

with(
  aa.cost,
  data.frame('Energy' = cost, row.names = AA)
) -> row.df

dat.mat |>
  pheatmap::pheatmap(
    scale = 'none',
    show_colnames = TRUE,
    show_rownames = TRUE,
    display_numbers = TRUE,
    cluster_cols = lib.clust,
    annotation_col = col.df,
    annotation_row = row.df,
    annotation_colors = cl,
    # fontsize = 1,
    number_color = 'black',
    color = colorRampPalette(rev(
      RColorBrewer::brewer.pal(n = 5, name = "RdBu")))(59),
    # filename = 'analysis/J_AA-expr-cor.jpeg', width = 8, height = 8
  ) -> heat.cor
dev.off()
  

################################################################################
# Produce 2 heatmaps, freqs and expression separatly but with the exact same
# row/gene order
################################################################################
# heatmap expression

# sneaky, add space to library name to give similar bottom height
sneak.len <-
  fz |>
  colnames() |>
  str_length() |>
  max()
sneak <- function(x) {
  # sprintf(paste0('%-', sneak.len + 6, 's'), x)
  x
}

dat.mat <- vz[my.mask, colnames(dat.cor)]
# rename to short sample name
look <-
  meta |>
  with(set_names(sneak(sample), lib))
colnames(dat.mat) <- look[colnames(dat.mat)]

cl <- list(
  'CO2' = meta |>
    pull(CO2) |>
    unique() |>
    sort() |>
    as.character() %>%
    set_names(cbPalette[c(1, 6, 2, 7)], .)
)
with(
  meta,
  data.frame('CO2' = as.character(CO2), row.names = sneak(sample))
) -> col.df


dat.mat |>
  pheatmap::pheatmap(
    scale = 'none',
    show_colnames = TRUE,
    show_rownames = FALSE,
    cutree_rows = 5,
    cluster_cols = heat.cor$tree_col,
    annotation_col = col.df,
    annotation_colors = cl,
    color = colorRampPalette(rev(
      RColorBrewer::brewer.pal(n = 5, name = "RdBu")))(59),
  ) -> heat.expr
dev.off()

################################################################################
# Heatmap AA freq

# apply(abs(dat.cor) >= .2, 1, any) |>
#   which() |>
#   names() -> aas
# dat.mat <- fz[my.mask, aas]

# dat.mat <- fz[my.mask, rownames(dat.cor)]
# 
# cl <- list(
#   'Energy' = c('white', 'black')
# )
# with(
#   aa.cost,
#   data.frame('Energy' = cost, row.names = AA)
# ) -> col.df
# 
# # rg <- max(abs(dat.mat)) is 8.6
# 
# dat.mat |>
#   pheatmap::pheatmap(
#     cluster_rows = heat.expr$tree_row,
#     cluster_cols = heat.cor$tree_row,
#     scale = 'none',
#     show_colnames = TRUE,
#     show_rownames = FALSE,
#     annotation_col = col.df,
#     annotation_colors = cl,
#     color = colorRampPalette(rev(
#       RColorBrewer::brewer.pal(n = 9 , name = "PuOr")))(39),
#       # RColorBrewer::brewer.pal(n = 11, name = "RdBu")))(39),
#     # breaks = seq(-rg, rg, length.out = 40),
#     # fix to keep extremes darker
#     breaks = c(
#       -8.7,
#       seq(-3, 3, length.out = 38),
#       8.7
#     )
#   ) -> heat.aa
# dev.off()

################################################################################

# !! Make a test run on expression the order really add up!
# dat.mat <- vz[my.mask, ]
dat.mat <- fz[my.mask, rownames(dat.cor)]

clusters <- cutree(heat.expr$tree_row, k = 5) |>
  as_tibble(rownames = 'Geneid') |>
  mutate(cluster = paste0('cluster', value))
    
clusters |>
  group_by(cluster) |>
  do(i = {
    dat.mat[.$Geneid, ] |>
      apply(2, mean) |>
      as_tibble(rownames = 'AA')
  }) |>
  ungroup() |>
  unnest(i) |>
  spread(AA, value) -> foo
dat.aa.clustered <-
  foo |>
  select(- cluster) |>
  as.matrix() |>
  magrittr::set_rownames(foo$cluster)

#!! enusre rows are in order as expected by tree
tibble(Geneid = with(heat.expr$tree_row, labels[order])) |>
  left_join(clusters, 'Geneid')  |>
  pull(cluster) |>
  unique() -> xs
dat.aa.clustered <- dat.aa.clustered[xs, ]

  
cl <- list(
  'Energy' = c('white', 'black')
)
with(
  aa.cost,
  data.frame('Energy' = cost, row.names = AA)
) -> col.df

rg <- max(abs(dat.aa.clustered))

dat.aa.clustered |>
  pheatmap::pheatmap(
    cluster_rows = FALSE,
    cluster_cols = heat.cor$tree_row,
    scale = 'none',
    show_colnames = TRUE,
    show_rownames = FALSE,
    annotation_col = col.df,
    annotation_colors = cl,
    color = colorRampPalette(rev(
      # RColorBrewer::brewer.pal(n = 9 , name = "PuOr")))(39),
      RColorBrewer::brewer.pal(n = 11, name = "RdBu")))(39),
    breaks = seq(-rg, rg, length.out = 40),
    # fix to keep extremes darker
    # breaks = c(
    #   -8.7,
    #   seq(-3, 3, length.out = 38),
    #   8.7
    # )
  ) -> heat.aa
dev.off()

################################################################################
  
cowplot::plot_grid(
  heat.cor$gtable,
  heat.expr$gtable,
  heat.aa$gtable,
  labels = 'AUTO',
  nrow = 1
)
ggsave(
  'analysis/J_AA-expr-cor.jpeg',
  width = 8 * 3, height = 8, dpi = 400
)


################################################################################
sessionInfo()
