# Investigate the impact of differential expression on amino acid compositions
# Part 1: Sanitize expression and establish independence of gene length
# Part 2: Transform AA frequencies for linear regression/Pearson
# Part 3: Check for correlations

library(tidyverse)
library(ggpubr)
library(patchwork)

set.seed(123)

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

expression <-
  'analysis/D_normalized-counts.tsv' |>
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

# gene lengths of coding genes
gene.lengths <-
  'raw-data/PCC7002-genome.gff.gz' |>
  rtracklayer::import.gff3() |>
  as_tibble() |>
  filter(gene_biotype == 'protein_coding') |>
  select(Geneid = ID, len = width)

deg30 <-
  'analysis/M_logFC-vs-30.tsv' |>
  read_tsv()

################################################################################
################################################################################
################################################################################
# Part 1: Sanitize expression and establish independence of gene length


expression |>
  pivot_longer(- Geneid) |>
  group_by(Geneid) |>
  summarize(avg = mean(value)) |>
  left_join(gene.lengths, 'Geneid') |>
  mutate_at('avg', ~ .x + 1) |>
  ggscatter(
    'avg', 'len',
    alpha = .6,
    add = 'reg.line',
    add.params = list(color = 'blue'),
    cor.coef = TRUE,
    cor.coeff.args = list(color = 'blue', size = 5)
  ) +
  scale_y_log10() +
  scale_x_log10() +
  xlab('Average expression, size-factor normalized') +
  ylab('Gene length, nt') +
  annotation_logticks() +
  theme_pubr(18) -> p
p.raw <- ggExtra::ggMarginal(p, fill = 'grey')

expression |>
  pivot_longer(- Geneid) |>
  group_by(Geneid) |>
  summarize(avg = mean(value)) |>
  left_join(gene.lengths, 'Geneid') |>
  mutate_at('avg', ~ .x + 1) |>
  mutate(ratio = avg/len) |>
  ggscatter(
    'ratio', 'len',
    alpha = .6,
    add = 'reg.line',
    add.params = list(color = 'blue'),
    cor.coef = TRUE,
    cor.coeff.args = list(color = 'blue', size = 5)
  ) +
  scale_y_log10() +
  scale_x_log10() +
  xlab('Average expression, relative to gene length') +
  ylab('Gene length, nt') +
  annotation_logticks() +
  theme_pubr(18) -> p
p.rel <- ggExtra::ggMarginal(p, fill = 'grey')


################################################################################
# Length norm expression matrix

rel.expression <-
  expression |>
  pivot_longer(- Geneid) |>
  inner_join(gene.lengths, 'Geneid') |>
  mutate(value = value / len) |>
  select(- len) |>
  drop_na() |>
  pivot_wider()

rel.mat <-
  rel.expression |>
  select(- Geneid) |>
  as.matrix() |>
  magrittr::set_rownames(rel.expression$Geneid)

# rel.mat <- vst.mat

################################################################################
# Justify cutoff for lowly expressed, not expression changing genes

# max logFC and avg rel expression
dat <-
  deg |>
  filter(is.de) |>
  group_by(Geneid) |>
  summarize(max.alf = max(abs(log2FoldChange))) |>
  ungroup() |>
  left_join(
    rel.expression |>
      pivot_longer(- Geneid) |>
      group_by(Geneid) |>
      summarize(avg = mean(value)),
    'Geneid'
  ) |>
  drop_na()

# top x% with some rounding on cutoff
dat |>
  select_if(is.numeric) |>
  # map(quantile, probs = seq(0, 1, .1))
  map2(c(.8, .8), ~ quantile(.x, probs = .y)) |>
  map(~ set_names(.x, NULL)) |>
  unlist() |>
  round() |>
  as.list() -> xs
# xs$baseMean <- round(xs$baseMean, -2) 
# xs$avg <- round(xs$avg, -2) 
xs

dat |>
  ggplot(aes(avg, max.alf)) +
  geom_point(alpha = .5) +
  geom_hline(yintercept = xs$max.alf, color = 'red') +
  geom_vline(xintercept = xs$avg, color = 'blue') +
  scale_y_continuous(breaks = 0:8) +
  scale_x_log10(labels = scales::label_comma()) +
  annotate('text', x = xs$avg - 1, y = 8, size = 8,
           hjust = 1,
           label = xs$avg, color = 'blue') +
  annotate('text', x = 500, y = 0.5, size = 8,
           label = xs$max.alf, color = 'red') +
  annotation_logticks(sides = 'b') +
  xlab('Average expression, relative to gene length') +
  ylab('Max. sig. abs. logFC') +
  theme_pubr(18)  -> p
  
p <- ggExtra::ggMarginal(p, fill = 'grey')

wrap_plots(p.raw) + wrap_plots(p.rel) + p +
  plot_annotation(tag_levels = 'A')
ggsave('analysis/J_expression-length-logFC.jpeg',
       width = 16, height = 7, dpi = 400)

################################################################################
# Expression mask

mask.deg <-
  dat |>
  filter(
    max.alf >= xs$max.alf,
    avg >= xs$avg
  ) |>
  pull(Geneid) |>
  unique()

length(mask.deg)
# 300

deg |>
  filter(is.de, abs(log2FoldChange) >= 1) |>
  pull(Geneid) %>%
  unique -> mask.deg

deg30 %>%
  filter(padj <= 0.001) |>
  filter(abs(log2FoldChange) >= 1) |>
  pull(Geneid) %>%
  unique -> mask.deg

################################################################################
################################################################################
################################################################################
# Part 2: Transform AA frequencies for linear regression/Pearson

# Check overal AA freq distributions

freqs.mat <-
  freqs |>
  select(- Geneid) |>
  as.matrix() |>
  magrittr::set_rownames(freqs$Geneid)

freqs |>
  pivot_longer(- Geneid) |>
  mutate(AA = fct_reorder(name, value)) -> aa.ordered
aa.ordered |>
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
# substantiate transformation

# robust z-scaler on a vector
myz <- function(x) {(x - median(x)) / mad(x) }

# scale per column/AA in a matrix
per.AA <- function(x.mat) {
  x.mat |>
    apply(2, myz) |>
    magrittr::set_rownames(rownames(x.mat))
}

# BoxCox on vector (matrix are converted to vectors)
vec.box <- function(x) {
  x <- c(x) - min(x) + 1
  # transform
  lmb <-forecast::BoxCox.lambda(x)
  forecast::BoxCox(x, lmb)
}

# BoxCox on matrix
my.box <- function(x.mat) {
  x2 <- vec.box(x.mat)
  # convert back to matrix
  x2 |>
    matrix(nrow = nrow(x.mat), ncol = ncol(x.mat)) |>
    magrittr::set_rownames(rownames(x.mat)) |>
    magrittr::set_colnames(colnames(x.mat))
}

# BoxCox per column/AA
box.AA <- function(x.mat) {
  x.mat |>
    apply(2, vec.box) |>
    magrittr::set_rownames(rownames(x.mat))
}
  
qs <- seq(0, 1, length.out = 1e5)
list(
  'Frequencies' = freqs.mat,
  'Z'           = freqs.mat             |> myz(),
  'aaZ'         = freqs.mat             |> per.AA(),
  'BC + Z'      = freqs.mat |> my.box() |> myz(),
  'BC + aaZ'    = freqs.mat |> my.box() |> per.AA(),
  # 'aaBC + Z'    = freqs.mat |> box.AA() |> myz(),
  # 'aaBC + aaZ'  = freqs.mat |> box.AA() |> per.AA(),
  'aaZ + BC + Z' = freqs.mat |> per.AA() |> my.box() |> myz()
) |>
  map(c) %>%
  map2(names(.), ~ tibble(x = qs, y = quantile(.x, probs = qs),
                          f = .y)) |>
  bind_rows() |>
  mutate_at('f', fct_inorder) |>
  mutate(x2 = qnorm(x)) -> dist.dat 

dist.dat |>
  filter(is.finite(x2)) |>
  mutate(error = y - x2) |>
  group_by(f) |>
  summarize(sse = sum(error ** 2)) |>
  arrange() |>
  arrange(sse) |>
  mutate(
    lab = sprintf("SSE: %.1f", sse),
    off = 1:n()
  ) -> sse.dat

dist.dat |>
  ggplot(aes(x2, y, color = f, group = f)) +
  geom_point() +
  geom_abline(slope = 1, linetype = 'dashed', color = 'black') +
  scale_color_manual(
    values = RColorBrewer::brewer.pal(9, 'Paired'),
    name = NULL
  ) +
  geom_text(
    aes(x = -2, y = 13 - off * 2, label= lab),
    size = 8,
    data = sse.dat, 
    show.legend = FALSE
  ) +
  xlab('Theoretic quantiles\nNormal distribution') +
  ylab('Observed quantiles')  +
  guides(color = guide_legend(ncol = 2)) +
  theme_pubr(18) -> p1.trans


################################################################################
# show after transform boxplots

# transform to a more "normal" distribution, best version from above
fz <-
  freqs.mat |>
  per.AA() |>
  my.box() |>
  myz()

# fz <- freqs.mat

fz |>
  as_tibble(rownames = 'Geneid') |>
  pivot_longer(- Geneid) |>
  mutate(AA = fct_relevel(name, levels(aa.ordered$AA))) |>
  ggplot(aes(AA, value, fill = AA)) +
  # geom_violin() +
  # geom_boxplot(fill = 'white', width = .2) +
  geom_boxplot() +
  xlab(NULL) +
  ylab('per AA z-scaled + BoxCox transformed\n& globally z-scaled') +
  theme_pubr(18) +
  theme(
    legend.position = 'hide',
    axis.text.x = element_text(angle = 60, hjust = 1)
  ) -> p1.z


  
(p1.total | (wrap_elements(full = p1.trans) + theme_pubr(18)) | p1.z) &
  plot_annotation(tag_levels = 'A')

ggsave('analysis/J_freqs-overall.jpeg', width = 18, height = 7, dpi = 400)

################################################################################
################################################################################
################################################################################
# Part 3: Check for correlations
#


# Z-Scaled rel expression per-gene

# rz <-
#   rel.mat |>
#   apply(1, scale) |>
#   t() |>
#   magrittr::set_colnames(colnames(rel.mat))

rz <-
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
# Heatmap of AA and expression correlation

# fz <-
#   freqs.mat |>
#   apply(2, rank)

mask.deg <- intersect(mask.deg, rownames(fz))

crossing(
  AA = colnames(fz),
  lib = colnames(rz)
) |>
  group_by(AA, lib) |>
  do(r2 = cor.test(fz[mask.deg, .$AA], rz[mask.deg, .$lib])$estimate) |>
  ungroup() |>
  unnest(r2) |>
  spread(lib, r2) -> dat

dat.cor <- dat |>
  select(- AA) |>
  as.matrix() |>
  magrittr::set_rownames(dat$AA)

# filter out low correlations
# dat.cor |>
#   abs() %>%
#   apply(1, max) -> foo
# dat.cor <- dat.cor[foo >= 0.2, ]

################################################################################
# Produce multiple heatmaps, freqs and expression separatly but with the shared
# row/gene hclust

# heatmap of correlations AA freq and Expression

dat.mat <- dat.cor

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
      pull(lib)
  ) |>
  as.hclust()

cl <- list(
  'CO2' = meta |>
    pull(CO2) |>
    unique() |>
    sort() |>
    as.character() %>%
    set_names(cbPalette[c(1, 6, 2, 7)], .),
  'Energy' = c('lightgrey', 'black')
)
with(
  meta,
  data.frame('CO2' = as.character(CO2), row.names = lib)
) -> col.df

with(
  aa.cost,
  data.frame('Energy' = cost, row.names = AA)
) -> row.df

rg <- dat.mat |> abs() |> max()
dat.mat |>
  pheatmap::pheatmap(
    scale = 'none',
    show_colnames = FALSE,
    show_rownames = TRUE,
    # display_numbers = TRUE,
    cluster_cols = lib.clust,
    annotation_col = col.df,
    annotation_row = row.df,
    annotation_colors = cl,
    number_color = 'black',
    color = colorRampPalette(rev(
      RColorBrewer::brewer.pal(n = 9, name = "PuOr")))(59),
    breaks = seq(-rg, rg, length.out = 60),
  ) -> heat.cor
dev.off()
  

################################################################################
################################################################################
# Additional heatmaps of expression (vst) clusters to amino acid concentrations

# build matrix
vst.mat <-
  vst |>
  select(- Geneid) |>
  as.matrix() |>
  magrittr::set_rownames(vst$Geneid)

deg |>
  filter(is.de) |>
  pull(Geneid) |>
  unique() -> de.mask

# dat.mat <- vst.mat[mask.deg, colnames(dat.cor)]
dat.mat <- vst.mat[de.mask, colnames(dat.cor)]

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
  data.frame('CO2' = as.character(CO2), row.names = lib)
) -> col.df

# Cluster columns ahead of pheatmap to rotate tree nicely
lib.clust <-
  dat.mat |>
  t() |>
  dist() |>
  hclust() |>
  ape::as.phylo() |>
  ape::rotateConstr(
    meta |>
      arrange(CO2, sample) |>
      pull(lib)
  ) |>
  as.hclust()

# no row clusters
CUT <- 6

# rg <- dat.mat |> abs() |> max()
dat.mat |>
  pheatmap::pheatmap(
    scale = 'row',
    clustering_method = 'ward.D2',
    # clustering_method = '',
    show_colnames = FALSE,
    show_rownames = FALSE,
    cutree_rows = CUT,
    cluster_cols = lib.clust,
    annotation_col = col.df,
    annotation_colors = cl,
    color = colorRampPalette(rev(
      RColorBrewer::brewer.pal(n = 9, name = "RdBu")))(59),
    # breaks = seq(-rg, rg, length.out = 60),
  ) -> heat.expr
dev.off()

#----------------------------------------------------------------------
# re-run pheatmap with cluster coloration

clusters <-
  cutree(heat.expr$tree_row, k = CUT) |>
  as_tibble(rownames = 'Geneid') %>%
  # combine with row order to have "nicer" order of cluster numbers in the heatmap
  `[`(heat.expr$tree_row$order, ) |>
  mutate(
    new.value = cumsum(lag(value, default = 0) != value)
  ) |>
  mutate(cluster = paste('cluster', new.value)) %>%
  with(data.frame(cluster = cluster, row.names = Geneid))

cl$cluster <-
  RColorBrewer::brewer.pal(CUT, 'Paired') |>
  # ggsci::pal_igv()(CUT) |>
  set_names(clusters$cluster |> unique())



dat.mat |>
  pheatmap::pheatmap(
    scale = 'row',
    clustering_method = 'ward.D2',
    show_colnames = FALSE,
    show_rownames = FALSE,
    cutree_rows = CUT,
    cluster_cols = lib.clust,
    annotation_col = col.df,
    annotation_row = clusters, # this is the new line compared to above
    annotation_colors = cl,
    color = colorRampPalette(rev(
      RColorBrewer::brewer.pal(n = 9, name = "RdBu")))(59),
    # breaks = seq(-rg, rg, length.out = 60),
  ) -> heat.expr
dev.off()
    
################################################################################
################################################################################
# distribution of AA freqs per cluster


# combine with AA freqs and clusters
fz[, rownames(dat.cor)] |>
  as_tibble(rownames = 'Geneid') |>
  pivot_longer(- Geneid, names_to = 'AA') -> foo
  # left_join(
foo |>
  inner_join(
    clusters  |>
      as_tibble(rownames = 'Geneid'),
    'Geneid'
  ) |>
  arrange(cluster)  -> dat
dat |>
  bind_rows(foo) |>
  mutate_at('cluster', replace_na, 'All genes') -> dat

dat |>
  ggplot(aes(cluster, value, fill = fct_inorder(cluster))) +
  geom_violin() +
  geom_boxplot(width = .5, fill = 'white') +
  xlab(NULL) +
  ylab('AA freuency') +
  scale_fill_brewer(palette = 'Paired', name = NULL) +
  # ggsci::scale_fill_igv(name = NULL) +
  geom_hline(
    aes(yintercept = avg),
    data = dat |>
      filter(cluster == 'All genes') |>
      group_by(AA) |>
      summarize(avg = median(value)),
    # color = 'black',
    # visually tether line to color of ref group
    color = RColorBrewer::brewer.pal(CUT + 1, 'Paired') |> last()
  ) +
  stat_compare_means(
    ref.group = 'All genes',
    # comparisons = cmps,
    hide.ns = TRUE,
    label = 'p.signif',
    method = 't.test',
    color = 'darkred'

  ) +
  facet_wrap(~ AA, scales = 'free_y') +
  # ylim(c(-8, 8)) +
  # facet_wrap(~ AA) +
  theme_pubr(18) +
  theme(axis.text.x = element_blank()) -> p

annotate_figure(
  p,
  bottom = text_grob(
    paste(
      'T-test vs all genes; P-value:',
      # 'ns  not significant;',
      '*  < 0.05;',
      '**  < 0.01;',
      '***  < 0.001;',
      '****  < 0.0001',
      sep = ' '
    ),
    color = 'darkred',
    hjust = 1, x = .9,
    face = "italic", size = 14,
  )
) -> p.boxes


################################################################################

################################################################################

cowplot::plot_grid(
  cowplot::plot_grid(
    heat.cor$gtable,
    heat.expr$gtable,
    nrow = 2,
    label_size = 18,
    labels = c("A", "B")
  ),
  p.boxes,
  labels = c(NA, "C"),
  rel_widths = c(1, 3),
  label_size = 18,
  nrow = 1
)
ggsave(
  'analysis/J_AA-expr-cor.jpeg',
  width = 20, height = 14, dpi = 400
)


################################################################################
sessionInfo()

################################################################################
# # explore connection with deg30 logFC and distribution of the AAs
# 
# xs <- seq(0, 1, .025)
# fz |>
#   c() |>
#   quantile(xs)
#   plot(x = xs)
#   # abs() |>
#   ecdf() -> fooo
# plot(fooo)
# 
# ?fooo()
# ?ecdf
# 
# 
# fz |>
#   abs() |>
#   apply(2, max) -> barr
# hist(barr)
# hist(fz)
# summary(c(fz))
# 
# 
# fz |>
#   abs() |>
#   apply(1, max) -> baz
# hist(baz)
# bar <- order(baz, decreasing = TRUE)[seq_len(500)]
# 
# # foo <- prcomp(fz[mask.deg, ])
# # foo <- prcomp(fz[baz >= 2, ])
# foo <- prcomp(fz[bar, ])
# percentVar <- foo$sdev^2/sum(foo$sdev^2)
# percentVar * 100
# 
# foo %>%
#   `[[`('x') |>
#   as_tibble(rownames = 'Geneid') |>
#   left_join(deg30, 'Geneid') |>
#   # ggplot(aes(PC1, PC2, size = abs(log2FoldChange), color = -log(padj))) +
#   ggplot(aes(PC1, PC2, color = abs(log2FoldChange), size = -log(padj))) +
#   # ggplot(aes(log2FoldChange, PC2, color = -log(padj))) +
#   geom_point(alpha = .7) +
#   # scale_color_gradient2(mid = 'grey', low = 'blue', high = 'red') +
#   scale_color_viridis_c() +
#   theme_dark()
# 
# fz[, rownames(dat.cor)] |>
#   as_tibble(rownames = 'Geneid') |>
#   pivot_longer(- Geneid, names_to = 'AA')  |>
#   left_join(deg30, 'Geneid') |>
#   filter(padj <= 0.05) |>
#   ggplot(aes(log2FoldChange, value, color = - log(padj))) +
#   geom_point(alpha = .5) +
#   scale_color_viridis_c() +
#   facet_wrap(~ AA)
# 
# pheatmap::pheatmap(fz[baz >= 2, ], show_rownames = FALSE,
#                    # scale = 'column',
#                    breaks = seq(-2, 2, length.out = 99),
#                    clustering_method = 'ward.D2')

################################################################################
# Overview

p1a <-
  freqs |>
  pivot_longer(- Geneid) |>
  ggdensity('value') +
  xlab('AA freq.') +
  theme_pubr(18)

p1b <-
  freqs |>
  pivot_longer(- Geneid) |>
  ggdensity('value', color = 'name') +
  xlab('AA freq.') +
  scale_color_manual(
    values = colorRampPalette(
      RColorBrewer::brewer.pal(12, 'Paired')
    )(20),
    name = NULL
  ) +
  theme_pubr(18)


(p1a | p1b)

################################################################################
# Estimate Beta parameters 

mask <- intersect(rownames(freqs.mat), mask.deg)
mat <- freqs.mat[mask, ] + .Machine$double.eps

qs <- seq(0, 3, length.out = 50)
# plot(qs, dweibull(qs, 1.5, .2))
plot(qs, dlnorm(qs, -2, .5))

freqs.paras <-
  mat |>
  # apply(2, \(x) MASS::fitdistr(unlist(x), 'lognormal')) |>
                               # start = list(meanlog = -2, sdlog = .5) ))|>
                               # lower = rep(.Machine$double.eps, 2))) |>
                               # method = 'L-BFGS-B')) |>
  # apply(2, \(x) MASS::fitdistr(unlist(x), 'weibull',
  #                              start = list(shape = 1.5, scale = .2) ))|>
                               # lower = rep(.Machine$double.eps, 2))) |>
                               # method = 'L-BFGS-B')) |>
  apply(2, \(x) MASS::fitdistr(unlist(x), 'beta',
                               start = list(shape1 = 3, shape2 = 50),
                               lower = rep(.Machine$double.eps, 2))) |>
                               # method = 'L-BFGS-B')) |>
  # apply(2, \(x) MASS::fitdistr(unlist(x), 'gamma',
  #                              start = list(shape = 1, rate = 1),
  #                              lower = rep(.Machine$double.eps, 2),
  #                              method = 'L-BFGS-B')) |>
  map('estimate') %>%
  map2(names(.), function(x, y) { x$name <- y ; x } ) |>
  bind_rows() |>
  mutate(
    # for lognormal
    # expected = exp(meanlog + sdlog ** 2 / 2),
    # var = (
    #   exp(meanlog ** 2) - 1
    # ) * exp(
    #   2 * meanlog * sdlog ** 2
    # )
    # for weibull
    # expected = scale *
    #   gamma(1 + 1 / shape),
    # var = scale ** 2 * (
    #   gamma(1 + 2 / shape) -
    #   gamma(1 + 1 / shape) ** 2
    # )
    # for beta
    expected = shape1 / (shape1 + shape2),
    var = shape1 * shape2 / (
      (shape1 + shape2) ** 2 +
      (shape1 + shape2 + 1)
    )
    # for gamma
    # expected = shape / rate,
    # var = shape / (rate ** 2)
  ) |>
  left_join(
    tibble(
      name = colnames(mat),
      obs.mean = apply(mat, 2, mean),
      obs.var = apply(mat, 2, var)
    ),
    'name'
  )


p2a <-
  freqs.paras |>
  ggscatter('expected', 'obs.mean',
            add = 'reg.line', add.params = list(color = 'blue'),
            cor.coef = TRUE, cor.coeff.args = list(color = 'red', size = 5)) +
  xlab('Theoretic') +
  ylab('Observed') +
  ggtitle('Expected Value') +
  theme_pubr(18)
p2b <-
  freqs.paras |>
  ggscatter('var', 'obs.var',
            add = 'reg.line', add.params = list(color = 'blue'),
            cor.coef = TRUE, cor.coeff.args = list(color = 'red', size = 5)) +
  xlab('Theoretic') +
  ylab('Observed') +
  ggtitle('Variance') +
  theme_pubr(18)

(p2a | p2b)

################################################################################

freqs.qs <-
  freqs.paras |>
  group_by(name) |>
  do(pf = {
    grp <- .
    partial(qbeta, shape1 = grp$shape1, shape2 = grp$shape2)
    # partial(qgamma, shape = grp$shape, rate = grp$rate)
  }) |>
  with(set_names(pf, name))


crossing(
  name = colnames(mat),
  qs = seq(0, 1, length.out = 1000)
) |>
  group_by_all() |>
  summarize(
    theo = freqs.qs[[name]](qs),
    empir = quantile(mat[, name], probs = qs)
  ) |>
  ungroup() -> qq.dat

qq.dat |>
  ggscatter(
    'theo', 'empir',
    color = 'name'
  )

################################################################################

################################################################################

freqs.ps <-
  freqs.paras |>
  select(name, shape1, shape2) |>
  nest_by(name) |>
  with(set_names(data, name)) |>
  map(~  partial(pbeta, shape1 = .x$shape1, shape2 = .x$shape2))

freqs.ecdf <-
  mat |>
  colnames() %>%
  set_names(.) |>
  map(~ mat[, .x]) |>
  map(ecdf)

# chi2 test for distribution equality
# 1. bin [0, 1] interval
qs <- seq(0, 1, length.out = 10)
qq.dat <-
  # 2. count observed frequency
  mat  |>
  as_tibble(rownames = 'gene') |>
  pivot_longer(- gene) |>
  mutate(v2 = cut(value, qs)) |>
  count(name, v2, name = 'obs') |>
  complete(name, v2) |>
  mutate_at('obs', replace_na, 0) |>
  # 3. estimate number of expected
  separate(v2, c('low', 'up'), sep = ',', remove = FALSE) |>
  mutate_at(c('up', 'low'), str_remove_all, '[^.0-9]') |>
  mutate_at(c('up', 'low'), as.numeric) |>
  # head()
  # group_by(name, v2, obs) |>
  rowwise() |>
  do(
    name = .$name,
    obs = .$obs,
    up = {
      grp <- .
      freqs.ps[[grp$name]](grp$up)
    },
    low = {
      grp <- .
      freqs.ps[[grp$name]](grp$low)
    }
  ) |>
  unnest(c(name, obs, up, low)) |>
  # mutate(exp = round((up - low) * nrow(mat)))  |>
  mutate(exp = (up - low) * nrow(mat))

freqs.x2 <-
  qq.dat |>
  filter((obs != 0) | (exp != 0)) |>
  mutate_at('obs', as.double) |>
  select(name, obs, exp) |>
  group_by(name) |>
  summarize(x2 = chisq.test(cbind(obs, exp))$p.value)

freqs.x2

crossing(
  name = colnames(mat),
  qs = seq(0, .5, length.out = 500)
) |>
  group_by_all() |>
  do(
    theo = freqs.ps[[first(.$name)]](.$qs),
    empir = freqs.ecdf[[first(.$name)]](.$qs)
  ) |>
  unnest(c(theo, empir)) |>
  pivot_longer(- c(name, qs), names_to = 'dist') |>
  mutate_at('dist', fct_recode,
            'Fitted density' = 'theo',
            'ECDF' = 'empir') |>
  ggplot(aes(qs, value, group = dist, color = dist)) +
  geom_line(size = 1.5) +
  geom_label(
    aes(x = 0.3, y = .2,
        group = NULL, color = NULL,
        label = sprintf(
          'shape1: %.1f\nshape2: %.1f\nX2-Test P: %.1e',
          shape1, shape2,
          x2
        )),
    data = freqs.x2 |> left_join(freqs.paras, 'name'),
    show.legend = FALSE
  ) +
  ggsci::scale_color_jama(name = NULL) +
  xlab('AA freq.') +
  ylab('Cumulative probability') +
  facet_wrap(~ name) +
  theme_bw(18)
  
################################################################################



################################################################################
# vst |>
expression |>
  pivot_longer(- Geneid) |>
  left_join(meta, c('name' = 'lib')) |>
  group_by(Geneid, CO2) |>
  # summarize(value = mean(value)) |>
  summarize(value = mean(value) |> round()) |>
  ungroup() |>
  mutate(CO2 = paste0('co2_', CO2)) |>
  pivot_wider(names_from = 'CO2') -> foo

foo |>
  select(-Geneid) |>
  as.data.frame() |>
  magrittr::set_rownames(foo$Geneid) -> vco2

vco2 <- vco2 + 1

# vco2 |>
#   apply(1, scale) |>
#   t() |>
#   magrittr::set_colnames(colnames(vco2)) |>
#   as.data.frame() -> vz2


deg |>
  filter(is.de) |>
  pull(Geneid) |>
  unique() |>
  intersect(rownames(freqs.mat)) -> mask

crossing(
  i = colnames(freqs.mat),
  # link = c("logit", "probit", "cloglog", "cauchit", "log", "loglog")
  # link = c("logit", "probit")
  # link = c("logit", "probit", "cloglog", "log", 'identity')
  link = 'log'
) |>
  # filter(i == 'alanine') |>
  group_by_all() |>
  do(r2 = {
    grp <- .
    print(grp)
    vco2[mask, ] |>
    # vz2[mask, ] |>
      mutate(AA = freqs.mat[mask, grp$i] + .Machine$double.eps) -> dat
    
    foo <-
      betareg::betareg(AA ~ .,
                       data = dat, link = grp$link)
                       # data = dat, link = make.link(grp$link))
    summary(foo)$pseudo.r.squared
  }) |>
  ungroup() |>
  unnest(r2) -> link.test
  
link.test |>
  spread(link, r2) |>
  select(- i) |>
  as.matrix() |>
  boxplot()


################################################################################

crossing(
  i = colnames(freqs.mat),
  link =  "probit"
) |>
  # filter(i == 'alanine') |>
  group_by_all() |>
  do(mod = {
    grp <- .
    vco2[mask, ] |>
      mutate(AA = freqs.mat[mask, grp$i] + .Machine$double.eps) -> dat
    
    betareg::betareg(AA ~ .,
                     data = dat, link = make.link(grp$link))
  })  -> mods

library(broom)
mods  |>
  group_by(i) |>
  do(mod = .$mod |> first() |> tidy() |> select(-component)) |>
  unnest(mod) -> quick.test

################################################################################

quick.test |>
  filter(term != '(phi)') |>
  filter(term != '(Intercept)') |>
  ggplot(aes(term, estimate, ymin = estimate - std.error, ymax = estimate + std.error,
             color = p.value)) +
  geom_errorbar() +
  geom_point() +
  scale_color_viridis_c(direction = -1) +
  facet_wrap(~ i) +
  theme_dark(12)
