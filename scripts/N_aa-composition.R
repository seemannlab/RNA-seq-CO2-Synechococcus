# Investigate the impact of differential expression on amino acid compositions
# Part 1: General overview
# Part 2: Sanitize Negative Binomial distribution
# Part 3: Identify "extreme" genes
# Part 4: Correlated to DEG

library(tidyverse)
library(ggpubr)
library(patchwork)

set.seed(123)

# https://riptutorial.com/r/example/28354/colorblind-friendly-palettes
cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442",
               "#0072B2", "#D55E00", "#CC79A7")

# nicer amino acid names
data(aaMap, package = 'Biobase')

# highlight 90% confidence interval in the various plots
# -> find upper/lower 5%
conf.alpha <- (1 - .9)/2

################################################################################
# Load amino acid composition, but as absolute counts

aas <- Biostrings::readAAStringSet('analysis/G_peptides.faa')
# The frequencies
freqs <- Biostrings::alphabetFrequency(aas)
rownames(freqs) <- names(aas)

# exclude all 0s columns
mask <- apply(freqs, 2, max) > 0
freqs2 <- freqs[, mask]
# nicer AA names
colnames(freqs2) <- with(aaMap, set_names(name, let.1))[colnames(freqs2)]


freqs <-
  freqs2 |>
  as_tibble(rownames = 'Geneid')

rm(freqs2)

################################################################################
# As matrix

freqs.mat <-
  freqs |>
  select(- Geneid) |>
  as.matrix() |>
  magrittr::set_rownames(freqs$Geneid)

################################################################################
# Get peptide lengths

freqs.len <-
  freqs.mat |>
  apply(1, sum) |>
  as.data.frame() |>
  set_names('length') |>
  as_tibble(rownames = 'Geneid')


################################################################################
################################################################################
# Load remaining data data

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

freqs.rel <-
  'analysis/G_frequencies.tsv' |>
  read_tsv()

deg30 <-
  'analysis/M_logFC-vs-30.tsv' |>
  read_tsv() |>
  mutate(is.de =  (padj <= 0.001) & (abs(log2FoldChange) >= 1) ) |>
  left_join(annot, 'Geneid')

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
################################################################################
# Part 1: General overview over the data

# Check overall AA freq distributions

freqs.rel |>
  pivot_longer(- Geneid) |>
  mutate(AA = fct_reorder(name, value)) -> aa.ordered
aa.ordered |>
  left_join(aa.cost, 'AA') |>
  ggplot(aes(AA, value, fill = cost.rank)) +
  geom_violin() +
  geom_boxplot(fill = 'white', width = .2) +
  scale_fill_viridis_c(name = 'Biosyntheis cost rank (Akashi et al.)') +
  xlab(NULL) +
  ylab('Gene length normalized\nAmino acid frequency') +
  theme_pubr(18) +
  theme(
    # legend.position = 'hide',
    axis.text.x = element_text(angle = 60, hjust = 1)
  )
ggsave('analysis/N_aa-overall.jpeg', width = 12, height = 8, dpi = 400)

################################################################################
    
# Show (lack of) correlation to gene length
freqs |>
  pivot_longer(- Geneid, names_to = 'AA') |>
  left_join(freqs.len, 'Geneid') |>
  mutate_at('value', ~ .x + 1) |>
  ggscatter(
    'length', 'value',
    alpha = .5,
    add = 'reg.line', add.params = list(color = 'red'),
    cor.coef = TRUE, cor.coeff.args = list(color = 'red', size = 5)
  ) +
  scale_x_log10(breaks = c(100, 1000)) +
  scale_y_log10() +
  annotation_logticks() +
  xlab('Peptide length') +
  ylab('No. amino acid occurence\n(with pseudo count)') +
  theme_pubr(18) +
  facet_wrap(~ AA, scales = 'free_y')

ggsave('analysis/N_correlation-aa-length.jpeg',
       width = 14, height = 8, dpi = 400)

################################################################################

# Check distribution of lengths
qq.normal.helper <- function(xs) {
  tibble(
    xs = xs[order(xs)],
    ps = ppoints(xs),
    expected = qnorm(ps)
  )
}

freqs.len |>
  pull(length) |>
  log10() |>
  qq.normal.helper() |>
  ggplot(aes(expected, xs)) +
  geom_point() +
  geom_smooth(method = 'lm') +
  stat_cor(size = 5, color = 'blue') +
  ylab('log10 gene length quantiles') +
  xlab('Normal distribution quantiles') +
  theme_pubr(18)


ggsave('analysis/N_length-log-normal.jpeg', width = 4, height = 4, dpi = 400)

################################################################################


################################################################################
################################################################################
################################################################################
# Part 2: Sanitize Negative Binomial distribution


# Fit Binomial model per Max Likelihood estimate (MLE)

freqs |>
  pivot_longer(- Geneid, names_to = 'AA') |>
  left_join(freqs.len, 'Geneid') -> foo
foo |>
  group_by(AA) |>
  summarize(
    # This is the MLE
    phat = sum(value) / sum(length)
  ) |>
  right_join(foo, 'AA') |>
  select(Geneid, AA, value, length, phat) |>
  mutate(
    expected = length * phat,
    low = qbinom(conf.alpha, length, phat),
    up = qbinom(1 - conf.alpha, length, phat)
  ) -> binomial.model

# determine confidence interval
binomial.conf <-
  binomial.model |>
  select(AA, phat, length) |>
  unique() |>
  rowwise() |>
  mutate(
    low = qbinom(conf.alpha, length, phat),
    up = qbinom(1 - conf.alpha, length, phat)
  )


################################################################################
# Plot Binomial model

binomial.conf |>
  pivot_longer(c('low', 'up')) |>
  mutate(l2 = ifelse(name == 'up', -length, length)) |>
  arrange(AA, name, l2) -> foo


binomial.model |>
  ggplot(aes(x = length, group = 1)) +
  geom_point(aes(y = value, color = 'observed'), size = .5, alpha = .3) +
  geom_polygon(aes(y = value, fill = '90% confidence'),
               alpha = .5, data = foo) +
  geom_line(aes(y = expected, color = 'expected'), size = 1.2) +
  scale_color_manual(values = cbPalette[c(6, 1)], name = NULL) +
  scale_fill_manual(values = cbPalette[c(7, 6, 1)], name = NULL) +
  scale_x_log10(breaks = c(100, 1000)) +
  annotation_logticks(sides = 'b') +
  xlab('Peptide length') +
  ylab('Absolute AA abundence') +
  facet_wrap(~ AA,scales = 'free') +
  theme_pubr(18)

ggsave('analysis/N_binomial-model.jpeg', width = 16, height = 9)


################################################################################

# Convert per Amino acid to function that has only dispersion
# open for optimization
aa.fns <-
  binomial.model %>%
  plyr::dlply(plyr::`.`(AA), function(tbl) {
    function(disp) {
      tbl |>
        group_by_all() |>
        reframe(p = dnbinom(value, mu = expected, size = 1 / disp)) |>
        pull(p) |>
        log() |>
        sum() |>
        # convert to min problem for optim
        prod(-1)
    }
  })


# Find optimal non-negative dispersion parameters
aa.disp <-
  aa.fns |>
  map(function(f) {
    optim(.5, f, method = 'L-BFGS-B',
          upper = .Machine$double.xmax,
          lower = .Machine$double.xmin)
  }) |>
  map('par') %>%
  tibble(AA = names(.), disp = .) |>
  unnest(disp)

################################################################################
# build NB model and confidence interval

nb.model <-
  binomial.model |>
  left_join(aa.disp, 'AA')

nb.conf <-
  nb.model |>
  select(AA, phat, length, disp, mu = expected) |>
  unique() |>
  rowwise() |>
  mutate(
    low = qnbinom(conf.alpha,    size = 1 / disp, mu = mu),
    up = qnbinom(1 - conf.alpha, size = 1 / disp, mu = mu)
  )

################################################################################
# Plot Negative-Binomial model

nb.conf |>
  pivot_longer(c('low', 'up')) |>
  mutate(l2 = ifelse(name == 'up', -length, length)) |>
  arrange(AA, name, l2) -> foo


nb.model |>
  ggplot(aes(x = length, group = 1)) +
  geom_point(aes(y = value, color = 'observed'), size = .5, alpha = .3) +
  geom_polygon(aes(y = value, fill = '90% confidence'),
               alpha = .5, data = foo) +
  geom_line(aes(y = expected, color = 'expected'), size = 1.2) +
  scale_color_manual(values = cbPalette[c(6, 1)], name = NULL) +
  scale_fill_manual(values = cbPalette[c(7, 6, 1)], name = NULL) +
  scale_x_log10(breaks = c(100, 1000)) +
  annotation_logticks(sides = 'b') +
  # Show dispersion estimates
  geom_text(
    aes(x = 100, y = y2, label = lab),
    data = aa.disp |>
      mutate(lab = sprintf('Disp.: %.4f', disp)) |>
      left_join(nb.conf |> group_by(AA) |> summarize(y2 = max(up) * .8), 'AA')
  ) +
  xlab('Peptide length') +
  ylab('Absolute AA abundence') +
  facet_wrap(~ AA,scales = 'free') +
  theme_pubr(18)

ggsave('analysis/N_negative-binomial-model.jpeg', width = 16, height = 9)


################################################################################
################################################################################
################################################################################
# Part 3: Identify "extreme" genes

# Compute P-Value of per chance observing No. AAs
aa.ps <-
  nb.model |>
  group_by_all() |>
  reframe(cum.p = pnbinom(value, mu = expected, size = 1 / disp)) |>
  # *2 for two-sided test
  mutate(extreme.p = ifelse(cum.p > .5, 1 - cum.p, cum.p) * 2) |>
  # fdr adjust
  group_by(AA) |>
  mutate(fdr = p.adjust(extreme.p, 'fdr')) |>
  ungroup()

aa.ps |>
  write_tsv('analysis/N_aa-analysis.tsv')

################################################################################

# Volcano-like plot

aa.ps |>
  ggplot(aes(value - expected, -log10(fdr), color = log10(length))) +
  geom_point() +
  scale_color_viridis_c() +
  xlab('Difference Observed - Expected AA abundance') +
  ylab('-log10(P-Value, FDR adjusted)') +
  geom_hline(yintercept = -log10(0.05), color = 'red') +
  geom_hline(yintercept = -log10(0.1), color = 'orange') +
  facet_wrap(~ AA, scales = 'free') +
  theme_pubr(18)

ggsave('analysis/N_volcano-like.jpeg', width = 12, height = 8)

################################################################################

short.list <-
  aa.ps |>
  filter(fdr <= .1) |>
  select(Geneid, AA, value, expected, length, phat, disp, fdr) |>
  left_join(annot, 'Geneid')


################################################################################
################################################################################
################################################################################
# Part 4: Correlated to DEG


list(
  # 'DE30 coding gene' = deg30 |>
  #   filter(is.de, type == 'protein_coding') |>
  #   pull(Geneid) |>
  #   unique(),
  'DE coding gene' = deg |>
    filter(is.de, type == 'protein_coding') |>
    pull(Geneid) |>
    unique(),
  '"Extreme" AA composition' = short.list$Geneid
) |>
  venn::venn(zcolor = 'style',  ilcs = 1.5, sncs = 1.5)

################################################################################

aa.enrich <-
  short.list |>
  drop_na(old_locus_tag) |>
  pull(old_locus_tag) |>
  clusterProfiler::enrichKEGG(
    organism = 'syp',
    pvalueCutoff = 0.05,
    pAdjustMethod = "fdr"
  ) |>
  as_tibble() |>
  transmute(
    KEGG = ID,
    Pathway = str_remove(Description, ' - Picosynechococcus sp. PCC 7002$'),
    GeneRatio, BgRatio,
    Pvalue = pvalue,
    FDR = p.adjust,
    genes = geneID
  ) |>
  # Extract enrichment factor
  mutate(
    path.size = BgRatio |>
      str_remove('/.*') |>
      as.integer(),
    genes.with.signal.in.path = GeneRatio |>
      str_remove('/.*') |>
      as.integer(),
    genes.with.signal = GeneRatio |>
      str_remove('.*/') |>
      as.integer(),
    prop.all.genes.in.pathway = BgRatio |>
      map(~ parse(text = .x)) |>
      map(eval) |>
      unlist(),
    expected.genes = genes.with.signal * prop.all.genes.in.pathway,
    enrichment = genes.with.signal.in.path / expected.genes
  )


################################################################################
aa.enrich |>
  mutate(
    y = sprintf('%s (%s)',  Pathway, KEGG) |>
      fct_reorder(enrichment),
    nice = sprintf(
      'Pathway size: %g\nGenes with "extreme" AA: %g\nEnrichment: %.1f, FDR: %.1e',
      path.size,
      genes.with.signal.in.path, 
      enrichment, FDR
    ),
    hj = ifelse(enrichment > 5, 1, 0),
    x2 = enrichment + ifelse(enrichment > 5, -1, 1) * 0.5,
  ) |>
  ggplot(aes(enrichment, y, color = - log10(FDR),
             size = path.size, label = nice)) +
  # geom_point(size = 8) +
  geom_point() +
  scale_size(range = c(5, 20), name = 'Pathway size') +
  geom_text(aes(x = x2, hjust = hj), color = 'black', size = 5) +
  scale_color_viridis_c() +
  guides(color = guide_colorbar(barwidth = unit(5, 'cm'))) +
  ylab(NULL) +
  xlab('Enrichment ratio observed over expected genes in pathway') +
  theme_pubr(18) +
  theme(
    panel.grid.major.y = element_line(linetype = 'dotted')
  )


################################################################################
paths <-
  'analysis/I_gene2pathway.tsv' |>
  read_tsv()

short.list |>
  semi_join(deg |> filter(is.de)) |>
  left_join(paths) |>
  View()

foo <-
  short.list |>
  semi_join(deg |> filter(is.de)) |>
  pull(Geneid)

short.list |>
  semi_join(deg |> filter(is.de)) |>
  left_join(paths) |>
  filter(Title == 'Photosynthesis') |>
  pull(Geneid) -> mask


pheatmap::pheatmap(vz.mat[foo, ], show_rownames = FALSE)
pheatmap::pheatmap(vz.mat[mask, ])

################################################################################
################################################################################
################################################################################




################################################################################
################################################################################

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
  # select(name, shape, rate) |>
  nest_by(name) |>
  with(set_names(data, name)) |>
  map(~  partial(pbeta, shape1 = .x$shape1, shape2 = .x$shape2))
  # map(~  partial(pgamma, shape = .x$shape, rate = .x$rate))

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
  qs = seq(0, 1, length.out = 500)
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
          # shape, rate,
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
