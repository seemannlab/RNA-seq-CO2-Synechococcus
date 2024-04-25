# Investigate the impact of differential expression on amino acid compositions
# Part 1: Sanitize expression and establish independence of gene length
# Part 2: Transform AA frequencies for linear regression/Pearson
# Part 3: Check for correlations

library(tidyverse)
library(ggpubr)
library(patchwork)

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

################################################################################
################################################################################
################################################################################
# Part 1: Sanitize expression and establish independence of gene length


expression |>
  pivot_longer(- Geneid) |>
  group_by(Geneid) |>
  summarize(avg = mean(value)) |>
  left_join(gene.lengths, 'Geneid') |>
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

pheatmap::pheatmap(fz, show_rownames = FALSE)

################################################################################
################################################################################
################################################################################
# Part 3: Check for correlations
#


# Z-Scaled rel expression per-gene

rz <-
  rel.mat |>
  apply(1, scale) |>
  t() |>
  magrittr::set_colnames(colnames(rel.mat))

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
dat.cor |>
  abs() %>%
  apply(1, max) -> foo
dat.cor <- dat.cor[foo >= 0.2, ]

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
    display_numbers = TRUE,
    cluster_cols = lib.clust,
    annotation_col = col.df,
    annotation_row = row.df,
    annotation_colors = cl,
    number_color = 'black',
    color = colorRampPalette(rev(
      RColorBrewer::brewer.pal(n = 9, name = "RdYlBu")))(59),
    breaks = seq(-rg, rg, length.out = 60),
  ) -> heat.cor
dev.off()
  

################################################################################
# heatmap expression


dat.mat <- rz[mask.deg, colnames(dat.cor)]

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
CUT <- 3

rg <- dat.mat |> abs() |> max()
dat.mat |>
  pheatmap::pheatmap(
    scale = 'none',
    show_colnames = FALSE,
    show_rownames = FALSE,
    cutree_rows = CUT,
    cluster_cols = lib.clust,
    annotation_col = col.df,
    annotation_colors = cl,
    color = colorRampPalette(rev(
      RColorBrewer::brewer.pal(n = 9, name = "RdYlBu")))(59),
    breaks = seq(-rg, rg, length.out = 60),
  ) -> heat.expr
dev.off()

#----------------------------------------------------------------------
# re-run pheatmap with cluster coloration

clusters <- cutree(heat.expr$tree_row, k = CUT) |>
  as_tibble(rownames = 'Geneid') |>
  mutate(cluster = paste('cluster', value)) %>%
  with(data.frame(cluster = cluster, row.names = Geneid))

cl$cluster <-
  ggsci::pal_igv()(CUT) |>
  set_names(clusters$cluster |> unique())



dat.mat |>
  pheatmap::pheatmap(
    scale = 'none',
    show_colnames = FALSE,
    show_rownames = FALSE,
    cutree_rows = CUT,
    cluster_cols = lib.clust,
    annotation_col = col.df,
    annotation_row = clusters, # this is the new line compared to above
    annotation_colors = cl,
    color = colorRampPalette(rev(
      RColorBrewer::brewer.pal(n = 9, name = "RdYlBu")))(59),
    breaks = seq(-rg, rg, length.out = 60),
  ) -> heat.expr
dev.off()
    
################################################################################

# average z-scaled len rel expression
azr <-
  rz |>
  as_tibble(rownames = 'Geneid') |>
  pivot_longer(- Geneid, names_to = 'lib') |>
  left_join(meta, 'lib') |>
  group_by(Geneid, CO2) |>
  summarize(avg.rel = mean(value))


# combine with AA freqs and clusters
fz[, rownames(dat.cor)] |>
  as_tibble(rownames = 'Geneid') |>
  pivot_longer(- Geneid, names_to = 'AA') |>
  # left_join(
  #   azr, 'Geneid',
  #   relationship = "many-to-many"
  # ) |>
  inner_join(
    clusters  |>
      as_tibble(rownames = 'Geneid'),
    'Geneid'
  ) -> dat

dat %>%
  pull(cluster) %>%
  unique %>%
  sort() %>%
  as.character() %>%
  crossing(a = ., b = .) |>
  filter(a < b) |>
  # tibble(a = ., b = lead(.)) %>%
  # drop_na() %>%
  rowwise() %>%
  do(i = unlist(.)) %>%
  pull(i) -> cmps

# dat$value |> summary()
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# -7.8724 -0.9708 -0.1180 -0.1368  0.7173  4.3109 

dat |>
  ggplot(aes(cluster, value, fill = cluster)) +
  geom_violin() +
  geom_boxplot(width = .5, fill = 'white') +
  xlab(NULL) +
  ylab('Tranformed AA freuency') +
  ggsci::scale_fill_igv() +
  stat_compare_means(
    comparisons = cmps,
    label = 'p.signif',
    method = 't.test',
    size = 5
    
  ) + 
  # facet_wrap(~ AA, scales = 'free_y') +
  ylim(c(-8, 8)) +
  facet_wrap(~ AA) +
  theme_pubr(18) +
  theme(axis.text.x = element_blank()) -> p

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
  labels = c(NA, "B"),
  label_size = 18,
  nrow = 1
)
ggsave(
  'analysis/J_AA-expr-cor.jpeg',
  width = 16, height = 12, dpi = 400
)


################################################################################
sessionInfo()

