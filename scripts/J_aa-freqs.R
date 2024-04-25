# Investigate the impact of differential expression on amino acid compositions
# Part 1: Sanitize expression and establish independence of gene length
# Part 2: Transform AA frequencies for linear regression/Pearson
# Part 3: Check for correlations

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
  round(1) |>
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
  annotate('text', x = 1e3, y = 0, size = 8,
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
# 173

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

# robust z-scaler
myz <- function(x) {(x - median(x)) / mad(x) }

qs <- seq(0, 1, length.out = 1e4)
list(
  'Frequencies' = freqs.mat,
  'Z-scaled' = freqs.mat |> myz(),
  'Sqrt' = sqrt(freqs.mat),
  'Sqrt & Z' = sqrt(freqs.mat) |> myz()
) |>
  map(c) %>%
  map2(names(.), ~ tibble(x = qs, y = quantile(.x, probs = qs),
                          f = .y)) |>
  bind_rows() |>
  mutate_at('f', fct_inorder) |>
  mutate(x2 = qnorm(x)) |>
  ggplot(aes(x2, y, color = f, group = f)) +
  geom_point() +
  geom_abline(slope = 1, linetype = 'dashed', color = 'black') +
  scale_color_manual(
    values = RColorBrewer::brewer.pal(4, 'Paired'),
    name = NULL
  ) +
  xlab('Theoreric quantile\nNormal distribution') +
  ylab('Observed quantile')  +
  theme_pubr(18) -> p1.trans

################################################################################
# transform to a more "normal" distribution


# fz <-
#   freqs.mat |>
#   apply(2, function(x) {(x - median(x)) / mad(x)}) |>
#   magrittr::set_rownames(rownames(freqs.mat))

fz <- freqs.mat |> sqrt() |> myz()

fz |>
  as_tibble(rownames = 'Geneid') |>
  pivot_longer(- Geneid) |>
  mutate(AA = fct_relevel(name, levels(aa.ordered$AA))) |>
  ggplot(aes(AA, value, fill = AA)) +
  # geom_violin() +
  # geom_boxplot(fill = 'white', width = .2) +
  geom_boxplot() +
  xlab(NULL) +
  ylab('Sqrt. Frequencies, Z-scaled') +
  theme_pubr(18) +
  theme(
    legend.position = 'hide',
    axis.text.x = element_text(angle = 60, hjust = 1)
  ) -> p1.z


(p1.total | p1.z | p1.trans) +
  plot_annotation(tag_levels = 'A')

ggsave('analysis/J_freqs-overall.jpeg', width = 17, height = 7, dpi = 400)

# pheatmap::pheatmap(fz, show_rownames = FALSE)


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
  'Energy' = c('lightgrey', 'black')
)
with(
  meta,
  data.frame('CO2' = as.character(CO2), row.names = sample)
) -> col.df

with(
  aa.cost,
  data.frame('Energy' = cost, row.names = AA)
) -> row.df

rg <- dat.mat |> abs() |> max()
dat.mat |>
  pheatmap::pheatmap(
    scale = 'none',
    show_colnames = TRUE,
    show_rownames = TRUE,
    display_numbers = TRUE,
    cluster_cols = lib.clust,
    cutree_rows = 3,
    annotation_col = col.df,
    annotation_row = row.df,
    annotation_colors = cl,
    number_color = 'black',
    color = colorRampPalette(rev(
      RColorBrewer::brewer.pal(n = 5, name = "RdYlBu")))(59),
    breaks = seq(-rg, rg, length.out = 60),
  ) -> heat.cor
dev.off()
  

################################################################################
# heatmap expression

# sneaky, add space to library name to give similar bottom height
sneak.len <-
  fz |>
  colnames() |>
  str_length() |>
  max()
sneak <- function(x) {
  sprintf(paste0('%-', sneak.len + 6, 's'), x)
  # x
}

dat.mat <- rz[mask.deg, colnames(dat.cor)]
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


rg <- dat.mat |> abs() |> max()
dat.mat |>
  pheatmap::pheatmap(
    scale = 'none',
    show_colnames = TRUE,
    show_rownames = FALSE,
    cutree_rows = 3,
    cluster_cols = heat.cor$tree_col,
    annotation_col = col.df,
    annotation_colors = cl,
    color = colorRampPalette(rev(
      RColorBrewer::brewer.pal(n = 5, name = "RdYlBu")))(59),
    breaks = seq(-rg, rg, length.out = 60),
  ) -> heat.expr
dev.off()

# TODO re-run pheatmap with cluster coloration

################################################################################
# Heatmap AA freq


dat.mat <- fz[mask.deg, rownames(dat.cor)]

cl <- list(
  'Energy' = c('lightgrey', 'black')
)
with(
  aa.cost,
  data.frame('Energy' = cost, row.names = AA)
) -> col.df


rg <- dat.mat |> abs() |> max()
dat.mat |>
  pheatmap::pheatmap(
    cluster_rows = heat.expr$tree_row,
    cluster_cols = heat.cor$tree_row,
    scale = 'column',
    cutree_rows = 3,
    cutree_cols = 3,
    show_colnames = TRUE,
    show_rownames = FALSE,
    annotation_col = col.df,
    annotation_colors = cl,
    color = colorRampPalette(rev(
      RColorBrewer::brewer.pal(n = 5, name = "RdYlBu")))(59),
    breaks = seq(-rg, rg, length.out = 60),
  ) -> heat.aa
dev.off()

################################################################################

# !! Make a test run on expression the order really add up!
# dat.mat <- rz[my.mask, ]
dat.mat <- fz[mask.deg, rownames(dat.cor)]

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
  'Energy' = c('lightgrey', 'black')
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
      RColorBrewer::brewer.pal(n = 9 , name = "PuOr")))(39),
      # RColorBrewer::brewer.pal(n = 9, name = "RdBu")))(39),
    # breaks = seq(-rg, rg, length.out = 40),
    # fix to keep extremes darker
    breaks = c(
      -rg,
      seq(-rg / 2, rg /2, length.out = 38),
     rg 
    )
  ) -> heat.aa2
dev.off()

################################################################################
  
cowplot::plot_grid(
  heat.cor$gtable,
  heat.expr$gtable,
  heat.aa$gtable,
  heat.aa2$gtable,
  labels = 'AUTO',
  nrow = 1
)
ggsave(
  'analysis/J_AA-expr-cor.jpeg',
  width = 8 * 3, height = 8, dpi = 400
)


################################################################################
sessionInfo()
