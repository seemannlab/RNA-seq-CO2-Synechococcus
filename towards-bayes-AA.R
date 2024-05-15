# Towards testing if AA composition influences the expression, test fit
# of various expression models against the RNA-seq data
# - Arithmatic average
# - DESeq expected expression

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

freqs <-
  'analysis/G_frequencies.tsv' |>
  read_tsv()


freqs.mat <-
  freqs |>
  select(- Geneid) |>
  as.matrix() |>
  magrittr::set_rownames(freqs$Geneid)

################################################################################
# General distribution overview

freqs |>
  pivot_longer(- Geneid) |>
  mutate(AA = fct_reorder(name, value)) -> aa.ordered

aa.colors <-
  aa.ordered |>
  pull(AA) |>
  levels() %>%
  set_names(
    colorRampPalette(RColorBrewer::brewer.pal(11, 'Spectral'))(length(.)),
    .
  )

aa.ordered |>
  ggplot(aes(AA, value, fill = AA)) +
  geom_violin() +
  geom_boxplot(fill = 'white', width = .2) +
  scale_fill_manual(values = aa.colors, name = NULL) +
  xlab(NULL) +
  ylab('Gene length normalized\nAmino acid frequency') +
  theme_pubr(18) +
  theme(
    legend.position = 'hide',
    axis.text.x = element_text(angle = 60, hjust = 1)
  ) -> p1.total

p1.total

################################################################################
################################################################################
################################################################################
# Part 1
# Question: Are frequencies beta distributed?
################################################################################

# given a vector of observations maximum likelihood estimate parameters
my.fitter <- function(xs) {
  # avoid absolute zero
  # double.xmin is "the smallest non-zero normalized floating-point number"
  xs2 <- xs + .Machine$double.xmin
  mle <- MASS::fitdistr(
    xs2, 'beta',
    start = list(shape1 = .3, shape2 = 4),
    lower = rep(.Machine$double.xmin, 2)
  )
  res <- mle$estimate |> as.list()
  # compute model expected variance/mean
  res$mean <- with(res, shape1 / (shape1 + shape2))
  res$var = with(res, shape1 * shape2 / (
    (shape1 + shape2) ** 2 +
    (shape1 + shape2 + 1)
  ))
  # compare with empirically observed variance/mean
  res$obs.mean = mean(xs)
  res$obs.var = var(xs)
  as_tibble(res)
}

# MLE estimate shapes per AA and globally
beta.fits <-
  freqs.mat |>
  apply(2, my.fitter) |>
  bind_rows() |>
  mutate(AA = colnames(freqs.mat))

################################################################################

# Plot overview of fitted/observed parameters
p1A <-
  beta.fits |>
  mutate_at('AA', fct_relevel, levels(aa.ordered$AA)) |>
  ggscatter('mean', 'obs.mean', color = 'AA',
            add = 'reg.line', add.params = list(color = 'blue'),
            cor.coef = TRUE, cor.coeff.args = list(color = 'red', size = 5)) +
  xlab('Fitted Beta distribution') +
  ylab('Observed') +
  scale_color_manual(values = aa.colors, name = NULL) +
  ggtitle('Expected Value') +
  theme_pubr(18)
p1B <-
  beta.fits |>
  mutate_at('AA', fct_relevel, levels(aa.ordered$AA)) |>
  ggscatter('var', 'obs.var', color = 'AA',
            add = 'reg.line', add.params = list(color = 'blue'),
            cor.coef = TRUE, cor.coeff.args = list(color = 'red', size = 5)) +
  scale_color_manual(values = aa.colors, name = NULL) +
  xlab('Fitted Beta distribution') +
  ylab('Observed') +
  ggtitle('Variance') +
  theme_pubr(18)

(p1.total /  (p1A | p1B)) / guide_area() +
  plot_annotation(tag_levels = 'A') +
  plot_layout(heights = c(3, 3, 1), guides = 'collect')

ggsave('~/Downloads/foo1.jpeg', width = 14, height = 12)

################################################################################

# xs <- seq(-2, 2, length.out = 100)
# plot(xs, dnorm(xs))
# plot(xs, pnorm(xs))
# plot(xs, qnorm(xs))

# Adapt from stats:::qqnorm.default to make a quantile-quantile plot
# for a fitted beta distribution

# For an observation vector and model parameters calc qqplot data
qq.beta.helper <- function(xs, shape1, shape2) {
  tibble(
    xs = xs[order(xs)],
    ps = ppoints(xs),
    expected = qbeta(ps, shape1, shape2)
  ) |>
    group_by(xs) |>
    slice_min(ps) |>
    ungroup()
}

qq.data <-
  freqs.mat |>
  colnames() |>
  map(function(i) {
    beta.fits |>
      filter(AA == i) |>
      with(qq.beta.helper(
        freqs.mat[, i],
        shape1,
        shape2
      )) |>
      mutate(AA = i)
  }) |>
  bind_rows()

p.qq1 <-
  qq.data |>
  mutate_at('AA', fct_relevel, levels(aa.ordered$AA)) |>
  ggplot(aes(expected, xs, color = AA)) +
  geom_point() +
  scale_color_manual(values = aa.colors, name = NULL) +
  geom_abline(slope = 1, linetype = 'dashed') +
  theme_pubr(18) +
  xlab('Fitted beta quantile') +
  ylab('Empirical quantile') +
  ggtitle('Quantile-quantile plot')

freqs.mat |>
  colnames() |>
  tibble(
    AA = _
  ) |>
  group_by(AA) |>
  reframe(
     aa.min = freqs.mat[, AA] |> min(),
     aa.max = freqs.mat[, AA] |> max()
  ) -> freqs.ranges

p.qq2 <-
  qq.data |>
  left_join(freqs.ranges, 'AA') |>
  filter(expected < aa.max) |>
  mutate_at('AA', fct_relevel, levels(aa.ordered$AA)) |>
  ggplot(aes(expected, xs, color = AA)) +
  geom_point() +
  scale_color_manual(values = aa.colors, name = NULL) +
  geom_abline(slope = 1, linetype = 'dashed') +
  theme_pubr(18) +
  xlab('Fitted beta quantile') +
  ylab('Empirical quantile') +
  ggtitle('Quantile-quantile plot',
          'Excluding fitted values above max. observed frequency')

((p.qq1 | p.qq2)) / guide_area() +
  plot_annotation(tag_levels = 'A') +
  plot_layout(heights = c(3, 1), guides = 'collect')

ggsave('~/Downloads/foo2.jpeg', width = 20, height = 7)

################################################################################
# Inspect half-normal plot of residuals


# For an observation vector and model parameters calc qqplot data
qq.normal.helper <- function(xs) {
  tibble(
    xs = xs[order(xs)],
    ps = ppoints(xs),
    expected = fdrtool::qhalfnorm(ps)
  )
}

# Compute residuals from beta model
residuals.helper <- function(xs, shape1, shape2) {
  mu <-  shape1 / (shape1 + shape2)
  abs(xs - mu)
}


qq.res.data <-
  freqs.mat |>
  colnames() |>
  map(function(i) {
    beta.fits |>
      filter(AA == i) |>
      with(residuals.helper(
        freqs.mat[, i],
        shape1,
        shape2
      )) |>
      qq.normal.helper() |>
      mutate(AA = i)
  }) |>
  bind_rows()

# simuate data for envelope
n <- nrow(freqs.mat)
nsim <- 100
# highlight 90% confidence interval
# -> find upper/lower 5%
alpha <- (1 - .9)/2
envelope.dat <-
  freqs.mat |>
  colnames() |>
  map(function(i) {
    # i <- 'asparagine'
    1:nsim |>
      # (it is ok if individual simulation runs fail)
      map(safely(function(j) {
        sampled <-
          beta.fits |>
          filter(AA == i) |>
          with(rbeta(n, shape1, shape2))
        sampled |>
          my.fitter() |>
          with(residuals.helper(sampled, shape1, shape2)) |>
          qq.normal.helper() |>
          mutate(AA = i, simulation = j)
      })) |>
      map('result') |>
      bind_rows()
  }) |>
  bind_rows()

# Estimate confidence intvervals
conf.band <-
  envelope.dat |>
  group_by(AA, expected) |>
  summarize(
    conf.up = quantile(xs, alpha),
    conf.low = quantile(xs, 1 - alpha),
    sim.runs = n()
  ) |>
  ungroup() |>
  pivot_longer(contains('conf')) |>
  left_join(freqs.ranges, 'AA') |>
  filter(expected < aa.max)



qq.res.data |>
  mutate_at('AA', fct_relevel, levels(aa.ordered$AA)) |>
  left_join(freqs.ranges, 'AA') |>
  filter(expected < aa.max) |>
  ggplot(aes(expected, xs, color = AA)) +
  geom_point() +
  geom_line(
    aes(y = value, group = paste(AA, name)),
    color = 'black',
    data = conf.band
  ) +
  scale_color_manual(values = aa.colors, name = NULL) +
  theme_pubr(18) +
  xlab('Half normal quantile') +
  ylab('Empirical quantile') +
  facet_wrap(~ AA) +
  theme(legend.position = 'hide') +
  ggtitle('Half-normal plot with simulated envelope')

ggsave('~/Downloads/foo2.jpeg', width = 12, height = 7)

################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
# Load data

annot <-
  'data/C_annotation.tsv' |>
  read_tsv()

meta <-
  'data/C_meta.tsv' |>
  read_tsv() |>
  mutate_at('CO2', ~ fct_reorder(as.character(.x), as.numeric(.x)))

raw.counts <-
  'data/C_raw-counts.tsv' |>
  read_tsv()

deg <-
  'analysis/D_stagewise-adjusted-DEGs.tsv' |>
  read_tsv()


################################################################################
# Build matrices


raw.counts.mat <-
  raw.counts %>%
  select(- Geneid) %>%
  as.matrix() %>%
  magrittr::set_rownames(raw.counts$Geneid)


################################################################################
# Overview

my.over <- function(mat) {
  tbl <-
    mat |>
    as_tibble(rownames = 'Geneid') |>
    pivot_longer(- Geneid)
  p1a <-
    tbl |>
    ggdensity('value') +
    xlab('AA freq.') +
    theme_pubr(18)
  
  p1b <-
    tbl |>
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
}

my.over(freqs.mat)
my.over(fz)

################################################################################

freqs.mat |>
  apply(2, function(x) {
    x <- x - min(x)
    # x <- x / max(x)
    x <- x / (max(x) + .Machine$double.eps * 10)
    # pmax(x, .Machine$double.xmin) |> pmin(1 - .Machine$double.eps * 10)
    x + .Machine$double.eps
  }) -> f2


my.over(f2)
boxplot(freqs.mat)
boxplot(f2)
summary(f2)

################################################################################

qs <- seq(0, 1, length.out = 50)
plot(qs, dbeta(qs, 5, 9))

# mat <- f2
mat <- freqs.mat + .Machine$double.xmin

freqs.paras <-
  mat |>
  apply(2, \(x) MASS::fitdistr(unlist(x), 'beta',
                               start = list(shape1 = 5, shape2 = 9),
                               lower = rep(.Machine$double.xmin, 2))) |>
  map('estimate') %>%
  map2(names(.), function(x, y) { x$name <- y ; x } ) |>
  bind_rows() |>
  mutate(
    # for beta
    expected = shape1 / (shape1 + shape2),
    var = shape1 * shape2 / (
      (shape1 + shape2) ** 2 +
      (shape1 + shape2 + 1)
    )
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

fz <-
  freqs.mat |>
  apply(2, scale) |>
  magrittr::set_rownames(rownames(freqs.mat))


################################################################################
# Exclude rRNA

annot %>%
  filter(type != 'rRNA') %>%
  pull(Geneid) -> mask

raw.noribo.mat <- raw.counts.mat[mask, ]

################################################################################
# Run DESeq2 to estimate the expected expression per condition

des <- DESeq2::DESeqDataSetFromMatrix(
  countData = raw.noribo.mat,
  colData = meta,
  design = ~ CO2 - 1
) |>
  DESeq2::DESeq()

################################################################################

# model below only for genes with expression changes
mask.de <-
  deg |>
  filter(is.de) |>
  select(Geneid) |>
  unique()

################################################################################
# Parameters of DESeq2 internal regression model

mod.des.expected <-
  des |>
  coef() |>
  as_tibble(rownames = 'Geneid') |>
  pivot_longer(- Geneid, names_to = 'CO2') |>
  mutate_at('CO2', str_remove, 'CO2')

mod.des.disps <-
  des |>
  DESeq2::dispersions() |> 
  set_names(rownames(des))

# Compute likelihood of NB model
# first, collect all parameters
mod.des.dat <-
  raw.noribo.mat |>
  as_tibble(rownames = 'Geneid') |>
  semi_join(mask.de, 'Geneid') |>
  pivot_longer(- Geneid, names_to = 'lib',  values_to = 'RNAseq') |>
  mutate(
    sf = DESeq2::sizeFactors(des)[lib],
    disp = mod.des.disps[Geneid]
  ) |>
  left_join(meta, 'lib') |>
  left_join(mod.des.expected, c('Geneid', 'CO2')) |>
  mutate(
    mu = sf * 2 ** value
  )

# sum up log probabilities for log likelihood
mod.des.ll <-
  mod.des.dat |>
  group_by_all() |>
  reframe(p = dnbinom(
    RNAseq,
    mu = mu,
    # for inverse of dispersion see DESeq2:::fitNbinomGLMs
    size = 1 / disp,
    log = TRUE
  )) |>
  pull(p) |>
  sum()

################################################################################

hist(freqs.mat[mask, 'alanine'])
hist(foo[mask, 'CO230'])

foo <- des |> coef()
bar <- des |> DESeq2::vst() |> SummarizedExperiment::assay()

qqnorm(foo[mask, ])
qqnorm(bar[mask, ])

crossing(
  aa = colnames(freqs.mat),
  co = colnames(foo)
) |>
  group_by_all() |>
  reframe(cor = cor.test(freqs.mat[mask, aa], foo[mask, co])$estimate) |>
  spread(co, cor) -> bar

bar |>
  select(- aa) |>
  as.matrix() |>
  magrittr::set_rownames(bar$aa)  |>
  pheatmap::pheatmap(display_numbers = TRUE)

################################################################################
# As a reference point simple linear model

des.norm <-  DESeq2::counts(des, normalized = TRUE)
des.norm <- log2(des.norm + 1)

example.lm <-
  lm(
    t(des.norm[mask.de$Geneid, ]) ~ CO2 - 1,
    data = meta
  )

# adapted from stat:::logLik.lm
# 0.5 * ( - N * (log(2 * pi) + 1 - log(N) +  log(sum(res^2))))
res <- example.lm$residuals |> c()
N <- length(res)
mod.lm.ll <-  0.5 * ( - N * (log(2 * pi) + 1 - log(N) +  log(sum(res^2))))

(- mod.des.ll) / mod.lm.ll

################################################################################
# arithmetic average per condition

mod.avg <-
  des.norm |>
  as_tibble(rownames = 'Geneid') |>
  pivot_longer(- Geneid, names_to = 'lib') |>
  left_join(meta, 'lib') |>
  group_by(Geneid, CO2) |>
  summarize(value = mean(value))

################################################################################
# quick scatter plot helper of observed expression vs various modeles


helper.join <- function(tbl, join) {
  des.norm |>
    as_tibble(rownames = 'Geneid') |>
    pivot_longer(- Geneid, names_to = 'lib',
                 values_to = 'obs') |>
    left_join(meta, 'lib') |>
    left_join(tbl, c('Geneid', join)) |>
    semi_join(mask.de, 'Geneid')
}
  
helper.scatter <- function(tbl, join, text) {
  tbl |>
    helper.join(join) |>
    ggscatter('value', 'obs',
              color = 'sample',
              add = 'reg.line', add.params = list(color = 'blue'),
              cor.coef = TRUE, cor.coeff.args = list(size = 5, color = 'red')) +
    ggsci::scale_color_igv(name = 'RNA-seq library') +
    ylab('Expression size-factor normalized, log2') +
    xlab(text) +
    theme_pubr(18)
}


################################################################################

helper.scatter(mod.avg, 'CO2', 'Per condition averaged expression') +
helper.scatter(mod.des.expected, 'CO2', 'DESeq2 expected expression')


bind_rows(
  mod.avg |>
    helper.join('CO2') |>
    mutate(
      diff = obs - value,
      mod = 'Per condition averaged expression'
    ),
  mod.des.expected |>
    helper.join('CO2') |>
    mutate(
      diff = obs - value,
      mod = 'DESeq2 expected expression'
    )
) |>
  ggecdf('diff', color = 'mod')

################################################################################
foo <-
  mod.des.expected |>
  semi_join(mask.de, 'Geneid') |>
  inner_join(
    freqs |>
      pivot_longer(- Geneid, names_to = 'AA', values_to = 'aa.freq'),
    'Geneid'
  )

foo |>
  ggplot(aes(value, aa.freq)) +
  geom_hex() +
  scale_fill_viridis_c() +
  stat_cor(color = 'white') +
  geom_smooth(method = 'lm', color = 'white', se = FALSE) +
  theme_dark() +
  facet_grid(AA ~ CO2, scales = 'free')

freqs |>
  pivot_longer(- Geneid, names_to = 'AA', values_to = 'aa.freq') |>
  # semi_join(mask.de, 'Geneid') |>
  mutate(value = mod.des.disps[Geneid]) |>
  ggplot(aes(log(value), aa.freq)) +
  geom_hex() +
  scale_fill_viridis_c() +
  stat_cor(color = 'white') +
  geom_smooth(method = 'lm', color = 'white', se = FALSE) +
  theme_dark() +
  facet_wrap(~ AA, scales = 'free')

################################################################################
################################################################################

mask <- intersect(
  mask.de$Geneid,
  rownames(freqs.mat)
)

Xs <- coef(des)[mask, ]
Xs <-
  Xs |>
  apply(1, scale) |>
  t() |>
  magrittr::set_colnames(colnames(Xs))
Xs <- Xs |> as.data.frame()

Ys <- freqs.mat[mask, ] + .Machine$double.xmin
Ys <- Ys |> as.data.frame()

# glm(Ys ~ ., data = Xs, family = 'Gamma')

foo <- lm(Xs ~ . - 1, data = Ys)


library(betareg)

AA.mod <-
  freqs.mat |>
  colnames() %>%
  set_names(., .) |>
  map(function(i) {
    Xs |>
      mutate(
        AA = freqs.mat[rownames(Xs), i],
        AA = AA + 1e-5
      ) |>
      summary()
      betareg(formula = AA ~ . - 1)
  })

################################################################################
################################################################################

library(broom)
AA.res <-
  AA.mod  |>
  map(tidy) %>%
  map2(names(.), ~ mutate(.x, AA = .y)) |>
  bind_rows()

################################################################################

foo |>
  tidy() |>
  # View()
# AA.res |>
  ggplot(aes(term, estimate,
             ymin = estimate - std.error, ymax = estimate + std.error,
             color = -log(p.value))) +
  geom_errorbar() +
  geom_point() +
  scale_color_viridis_c(direction = +1) +
  coord_flip() +
  facet_wrap(~ response) +
  # facet_grid(AA ~ component) +
  theme_dark(12)

################################################################################
################################################################################

des |>
  coef()  |>
# Xs |>
  as_tibble(rownames = 'Geneid') |>
  left_join(freqs, 'Geneid') |>
  ggplot(aes(CO28, lysine, color = lysine)) +
  scale_color_viridis_c() +
  geom_point() +
  theme_dark()


################################################################################
################################################################################


# foo <- prcomp(fz[bar, ])
foo <- prcomp(Xs[mask, ])
percentVar <- foo$sdev^2/sum(foo$sdev^2)
percentVar * 100

foo %>%
  `[[`('x') |>
  as_tibble(rownames = 'Geneid') |>
  left_join(freqs, 'Geneid') |>
  ggplot(aes(PC1, PC2, color = valine)) +
  geom_point(alpha = .7) +
  # scale_color_gradient2(mid = 'grey', low = 'blue', high = 'red') +
  scale_color_viridis_c() +
  theme_dark()

pheatmap::pheatmap(fz[baz >= 2, ], show_rownames = FALSE,
                   # scale = 'column',
                   breaks = seq(-2, 2, length.out = 99),
                   clustering_method = 'ward.D2')

################################################################################
################################################################################

foo <- 
  X |>
  mutate(aa = predict(AA.mod$alanine, X)) |>
  rename(expr = CO28)


tibble(
  aa = Ys[mask, 'valine'],
  expr = Xs[mask, 'CO230']
) |>
  ggplot(aes(expr, aa)) +
  geom_point(alpha = .2) +
  geom_density_2d(size = 1.2) +
  # geom_point(data = foo, color = 'red') +
  stat_cor(size = 5, color = 'red') +
  geom_smooth(method = 'lm', se = FALSE, color = 'red')


################################################################################
################################################################################

freqs |>
  pivot_longer(- Geneid, names_to = 'AA', values_to = 'aa.freq') |>
  inner_join(deg, 'Geneid') |>
  # filter(is.de) |>
  ggplot(aes(-log(padj), aa.freq)) +
  geom_hex() +
  scale_fill_viridis_c() +
  stat_cor(color = 'white') +
  geom_smooth(method = 'lm', color = 'white', se = FALSE) +
  theme_dark() +
  facet_grid(test ~ AA, scales = 'free')

################################################################################

mod.des.expected |>
  mutate(disp = mod.des.disps[Geneid]) |>
  ggscatter('value', 'disp', facet.by = 'CO2')

################################################################################

mod.des.expected

################################################################################
################################################################################
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


################################################################################
################################################################################
# Overview, after z-scale

p1a <-
  fz |>
  as_tibble(rownames = 'Geneid') |>
  pivot_longer(- Geneid) |>
  ggdensity('value') +
  xlab('AA freq.') +
  theme_pubr(18)

p1b <-
  fz |>
  as_tibble(rownames = 'Geneid') |>
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
# mat <- freqs.mat[mask, ] + .Machine$double.eps
mat <- freqs.mat + .Machine$double.eps

qs <- seq(0, 1, length.out = 50)
plot(qs, dbeta(qs, 3, 10))

freqs.paras <-
  mat |>
  apply(2, \(x) MASS::fitdistr(unlist(x), 'beta',
                               start = list(shape1 = 3, shape2 = 50),
                               lower = rep(.Machine$double.eps, 2))) |>
                               # method = 'L-BFGS-B')) |>
  map('estimate') %>%
  map2(names(.), function(x, y) { x$name <- y ; x } ) |>
  bind_rows() |>
  mutate(
    # for beta
    expected = shape1 / (shape1 + shape2),
    var = shape1 * shape2 / (
      (shape1 + shape2) ** 2 +
      (shape1 + shape2 + 1)
    )
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

mask <- intersect(
  rownames(freqs.mat),
  rownames(vst.mat)
)


vst.mat <-
  vst |>
  select(- Geneid) |>
  as.matrix() |>
  magrittr::set_rownames(vst$Geneid)

rz <-
  vst.mat |>
  apply(1, scale) |>
  t() |>
  magrittr::set_colnames(colnames(vst.mat))


i <- 'valine'
rz[mask, ] |>
  as_tibble(rownames = 'Geneid') |>
  pivot_longer(- Geneid, names_to = 'lib') |>
  left_join(meta) |>
  mutate(y = freqs.mat[Geneid, i] + .Machine$double.eps) -> dat

# head(dat)
# Geneid                  lib                                      value CO2   sample lane       y
# gene-SYNPCC7002_RS00010 CO2-0-04percent_Lane1-S1_CO2-adaptation -0.591 0.04  S1     Lane1 0.0491
# gene-SYNPCC7002_RS00010 CO2-0-04percent_Lane2-S2_CO2-adaptation  0.611 0.04  S2     Lane2 0.0491
# gene-SYNPCC7002_RS00010 CO2-0-04percent_Lane2-S3_CO2-adaptation  1.51  0.04  S3     Lane2 0.0491
# gene-SYNPCC7002_RS00010 CO2-0-04percent_Lane2-S4_CO2-adaptation  1.04  0.04  S4     Lane2 0.0491
# gene-SYNPCC7002_RS00010 CO2-30percent_Lane1-S14_CO2-adaptation   0.384 30    S14    Lane1 0.0491
# gene-SYNPCC7002_RS00010 CO2-30percent_Lane1-S16_CO2-adaptation  -0.889 30    S16    Lane1 0.0491


my.mod <-
  betareg::betareg(
    y ~ CO2 : value | CO2,
    data = dat
  )

dat |>
  mutate(pred = predict(my.mod, dat)) |>
  ggplot(aes(value, y)) +
  geom_point(aes(color = lib), alpha = 0.5) +
  geom_density_2d(show.legend = FALSE) +
  geom_smooth(method = 'lm', se = FALSE, show.legend = FALSE) +
  # geom_point(aes(y = pred), color = 'red') +
  geom_point(aes(x = pred), color = 'red') +
  facet_wrap(~ CO2)


summary(my.mod)

lmtest::coeftest(my.mod)

glm(
  value ~ CO2 : y,
  family = 'gaussian',
  data = dat
) -> my.mod

################################################################################
# convert freqs to quantiles

fps <-
  freqs.paras |>
  group_by(name) |>
  do(pf = {
    grp <- .
    partial(pbeta, shape1 = grp$shape1, shape2 = grp$shape2)
  }) |>
  with(set_names(pf, name))


freqs |>
  pivot_longer(- Geneid) |>
  group_by_all() |>
  reframe(ps = fps[[name]](value)) -> foo
  # ggplot(aes(value, ps)) +
  # geom_point() +
  # facet_wrap(~ name)

foo |>
  ggplot(aes(ps, color = name)) +
  geom_density()

str(fps)

################################################################################



################################################################################

mat


################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################