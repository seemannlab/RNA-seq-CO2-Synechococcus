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


signals <-
  'analysis/I_signals.tsv' |>
  read_tsv()

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
  mutate(ratio = log2( (value + 1) / (expected + 1) )) |>
  ggplot(aes(ratio, -log10(fdr), color = log10(length))) +
  geom_point() +
  scale_color_viridis_c(name = 'log10 Peptide length') +
  xlab('log2 ratio Observed / Expected AA abundance') +
  ylab('-log10(P-Value, FDR adjusted)') +
  geom_hline(yintercept = -log10(0.05), color = 'red') +
  facet_wrap(~ AA, scales = 'free') +
  theme_pubr(18)

ggsave('analysis/N_volcano-like.jpeg', width = 12, height = 8)

################################################################################

short.list <-
  aa.ps |>
  filter(fdr <= .05) |>
  select(Geneid, AA, value, expected, length, phat, disp, fdr) |>
  left_join(annot, 'Geneid')


################################################################################
################################################################################
################################################################################
# Part 4: Correlated to DEG


list(
  'DE coding gene' = deg |>
    filter(is.de, type == 'protein_coding') |>
    pull(Geneid) |>
    unique(),
  '"Extreme" AA composition' = short.list$Geneid,
  'Secretion Signal' = signals |>
    filter(signal != 'No SP') |>
    pull(Geneid)
) |>
  venn::venn(zcolor = 'style',  ilcs = 1.5, sncs = 1.5,
             box = FALSE, ggplot = TRUE)

ggsave(
  'analysis/N_extreme-AA-deg.jpeg',
  width = 10, height = 10, dpi = 400
)

################################################################################
# Heatmap of expression

my.list <-
  short.list |>
  semi_join(deg |> filter(is.de), 'Geneid') |>
  transmute(
    Geneid,
    old_locus_tag, AA, observed = value, expected,
    ratio = log2( (value + 1) / (expected + 1) ),
    length, fdr,
    name, product
  ) |>
  arrange(desc(ratio)) |>
  mutate(
    old_locus_tag = ifelse(
      is.na(old_locus_tag),
      Geneid,
      old_locus_tag
    ) |> str_remove('^.*SYNPCC7002_'),
    nice = sprintf('%s (%s)',
                   str_replace(product, '(?<=.{45}).....+$', '...'),
                   old_locus_tag)
  )

# build matrix
vst.mat <-
  vst |>
  select(- Geneid) |>
  as.matrix() |>
  magrittr::set_rownames(vst$Geneid)

dat.mat <- my.list |>
  select(Geneid, old_locus_tag) |>
  unique() |>
  with(vst.mat[Geneid, ] |> magrittr::set_rownames(old_locus_tag))

with(
  meta,
  data.frame('CO2' = as.character(CO2), row.names = lib)
) -> col.df


col.row <-
  signals |>
  filter(signal != 'No SP') |>
  inner_join(my.list, 'Geneid') |>
  select(ol = old_locus_tag.y, signal) |>
  unique() |>
  with(data.frame(Signal = signal, row.names = ol))

cl <- list(
  'CO2' = meta |>
    pull(CO2) |>
    unique() |>
    sort() |>
    as.character() %>%
    set_names(cbPalette[c(1, 6, 2, 7)], .),
  Signal = col.row |>
    pull(Signal) |>
    unique() %>%
    set_names(
      RColorBrewer::brewer.pal(length(.), 'Dark2'),
      .
    )
)

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


dat.mat |>
  pheatmap::pheatmap(
    scale = 'row',
    show_colnames = FALSE,
    border_color = NA,
    cluster_cols = lib.clust,
    annotation_col = col.df,
    annotation_row = col.row,
    annotation_names_row = FALSE,
    annotation_colors = cl,
    color = colorRampPalette(rev(
      RColorBrewer::brewer.pal(n = 9, name = "RdBu")))(59),
  ) -> heat.expr
dev.off()

    
################################################################################
# AA ratios per gene and AA as colored raster/table

p.ratio <-
  my.list |>
  mutate_at('old_locus_tag', fct_relevel,
            with(heat.expr$tree_row, labels[order]) ) |>
  arrange(desc(old_locus_tag)) |>
  mutate_at('nice', fct_inorder) |>
  mutate(
    lab = round(ratio, 2),
    lab.cl = ifelse(lab < 0, 'white', 'black')
  ) |>
  ggplot(aes(AA, nice, fill = ratio, label = lab, color = I(lab.cl))) +
  geom_raster() +
  geom_text() +
  scale_fill_viridis_c(name = 'log2 ratio Observed / Expected AA abundance') +
  guides(fill = guide_colorbar(barwidth = unit(3, 'cm'))) +
  xlab(NULL) +
  ylab(NULL) +
  theme_pubr() +
  theme(axis.text.x = element_text(angle = 20, hjust = 1))

################################################################################

cowplot::plot_grid(
  heat.expr$gtable,
  p.ratio,
  labels = 'AUTO',
  rel_widths = c(1, 2),
  label_size = 18,
  nrow = 1
)
ggsave(
  'analysis/N_extreme-AA.jpeg',
  width = 14, height = 7, dpi = 400
)

################################################################################
################################################################################
sessionInfo()
