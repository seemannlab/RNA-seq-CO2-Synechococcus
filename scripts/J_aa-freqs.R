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
  geom_boxplot() +
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
  geom_hline(yintercept = c(-2, 2), linetype = 'dashed') +
  theme(
    legend.position = 'hide',
    axis.text.x = element_text(angle = 60, hjust = 1)
  ) -> p1.z

(p1.total | p1.z) +
  plot_annotation(tag_levels = 'A')

ggsave('analysis/J_freqs-overall.jpeg', width = 14, height = 7, dpi = 400)

################################################################################
# Revisit gene length norm

'raw-data/PCC7002-genome.gff.gz' |>
  rtracklayer::import.gff3() |>
  as_tibble() |>
  filter(type == 'gene') |>
  select(Geneid = ID, len = width) -> lens

norm.counts <-
  'analysis/D_normalized-counts.tsv' |>
  read_tsv()

norm.counts |>
  pivot_longer(- Geneid, names_to = 'lib') |>
  left_join(meta, 'lib') |>
  group_by(Geneid, CO2) |>
  summarize(av = mean(value)) |>
  ungroup() |>
  inner_join(lens, 'Geneid') |>
  mutate(
    rel.av = av / len,
    rel.log.av = log1p(rel.av)
  ) -> rel.counts

vst.mat |>
  apply(1, mean) |>
  enframe(name = 'Geneid') |>
  left_join(lens) |>
  ggscatter('value', 'len', add = 'reg.line', cor.coef = TRUE)

################################################################################

mask2 <-
  fz |>
  abs() |>
  apply( 1, max) %>%
  `>=`(2) |>
#   table()
# FALSE  TRUE 
# 1731  1455 
  which() |>
  names()

vst |>
  pivot_longer(- Geneid, names_to = 'lib') |>
  left_join(meta, 'lib') |>
  group_by(Geneid, CO2) |>
  summarize(av = mean(value)) -> d1


deg |>
  filter(is.de, abs(log2FoldChange) >= 1) |>
  select(Geneid) |>
  unique() |>
  left_join(d1, 'Geneid') |>
  # left_join(rel.counts, 'Geneid') |>
  inner_join(
    # freqs |>
    fz |>
      as_tibble(rownames = 'Geneid') |>
      pivot_longer(- Geneid, names_to = 'AA', values_to = 'AA.freq'),
    'Geneid',
    relationship = "many-to-many"
  ) |>
  drop_na() -> d2

d2 |>
  filter(Geneid %in% mask) |>
  ggplot(aes(av, AA.freq)) +
  geom_hex(bins = 20) +
  scale_fill_viridis_c() +
  geom_smooth(method = 'lm', se = FALSE, color = 'red') +
  stat_cor(size = 5, color = 'red') +
  # ylim(c(-5, 5)) +
  # facet_wrap(~ AA)
  # facet_grid(AA ~ CO2, scales = 'free_y')
  facet_grid(CO2 ~ AA)

################################################################################

avg.vst <-
  vst |>
  pivot_longer(- Geneid, names_to = 'lib') |>
  left_join(meta, 'lib') |>
  group_by(Geneid, CO2) |>
  summarize(av = mean(value)) |>
  ungroup()

xcut <-
  avg.vst |>
  pull(av) |>
  quantile(.9) |>
  round(1)

avg.vst |>
  ggplot(aes(av, color = CO2)) +
  stat_ecdf(size = 2, alpha = .7) +
  scale_color_manual(values = cbPalette[c(1, 6, 2, 7)]) +
  xlab('Average gene expression (vst)') +
  ylab('ECDF') +
  geom_hline(yintercept = .9, color = 'red') +
  geom_vline(xintercept = xcut, color = 'blue') +
  annotate(
    'text', x = xcut + 2, y = .1,
    color = 'blue', size = 8,
    label = xcut
  ) +
  scale_y_continuous(breaks = seq(0, 1, .1)) +
  theme_pubr(18)

################################################################################


avg.vst |>
  # filter(av >= xcut) |>
  semi_join(
    deg |>
      filter(is.de) |>
      # filter(is.de, abs(log2FoldChange) >= 1) |>
      select(Geneid) |>
      unique(),
    'Geneid'
  ) |>
  left_join(
    fz |>
      as_tibble(rownames = 'Geneid') |>
      pivot_longer(- Geneid, names_to = 'AA'),
    'Geneid',
    relationship = "many-to-many"
  ) -> foo

foo |>
  drop_na() |>
  ggplot(aes(av, value)) +
  geom_hex() +
  geom_smooth(method = 'lm', color = 'blue') +
  stat_cor(color = 'red', size = 5) +
  scale_fill_viridis_c() +
  facet_grid(CO2 ~ AA)

foo |>
  drop_na() |>
  ggplot(aes(CO2, value, group = CO2)) +
  geom_violin() +
  facet_wrap(~ AA)

################################################################################

################################################################################
# Towards highly expressed genes

deg |>
  filter(is.de) |>
  group_by(Geneid, name, product, baseMean) |>
  summarize(max.alf = max(abs(log2FoldChange))) |>
  ungroup() |>
  # select_if(is.numeric) |>
  # map(quantile, probs = seq(0, 1, .1))
  # baseMean
  # 0%          10%          20%          30%          40%          50%          60%          70%          80%          90%         100% 
  # 1.016349e+01 2.017723e+02 3.723039e+02 5.683223e+02 8.075876e+02 1.112197e+03 1.558964e+03 2.208219e+03 3.682647e+03 8.170388e+03 2.175254e+06 
  # 
  # $max.alf
  # 0%       10%       20%       30%       40%       50%       60%       70%       80%       90%      100% 
  # 0.1457555 0.3764709 0.5005927 0.6157282 0.7303800 0.8547823 0.9970898 1.1576052 1.4102985 1.8326079 8.6502626 
  # -> the red lines below are approx the 20%
  ggplot(aes(baseMean, max.alf)) +
  # geom_hex(bins = 100) +
  # scale_fill_viridis_c() +
  geom_point(alpha = .5) +
  geom_hline(yintercept = 1, color = 'red') +
  geom_vline(xintercept = 400, color = 'blue') +
  scale_y_continuous(breaks = 0:8) +
  scale_x_log10(labels = scales::label_comma()) +
  annotate('text', x = 1000, y = 0, size = 8,
           label = '400', color = 'blue') +
  annotate('text', x = 10, y = 0, size = 8,
           label = '1', color = 'red') +
  xlab('Average expression') +
  ylab('Max. abs. logFC') +
  theme_pubr(18) -> p
  
p <- ggExtra::ggMarginal(p, fill = 'grey')
ggsave('analysis/J_expression-logFC.jpeg',plot = p,
       width = 9, height = 7, dpi = 400)

# Focus on DEGs
deg |>
  filter(is.de, abs(log2FoldChange) >= 1, baseMean >= 400) |>
  # filter(is.de) |>
  pull(Geneid) |>
  unique() |>
  intersect(rownames(fz)) -> mask

################################################################################

deg |>
  select(Geneid, test, log2FoldChange) |>
  spread(test, log2FoldChange) |>
  drop_na() |>
  filter(Geneid %in% mask) -> foo
foo |>
  select(-Geneid) |>
  as.matrix() |>
  magrittr::set_rownames(foo$Geneid) -> lfc.mat
################################################################################

# BiocManager::install('mixOmics')
# library(mixOmics)

rel.counts |>
  select(Geneid, CO2, rel.log.av) |>
  spread(CO2, rel.log.av) -> foo

foo |>
  select(- Geneid) |>
  as.matrix() |>
  magrittr::set_rownames(foo$Geneid) -> rel.mat

# Per gene z-Scaled expression
rz <-
  rel.mat |>
  apply(1, scale) |>
  t() |>
  magrittr::set_colnames(colnames(rel.mat))

mm <- intersect(mask, mask2)

# X <- rel.mat[mask, ]
X <- rz[mm, ]
Y <- fz[mm, ]
# X <- vst.mat[mask, ]

pls.result <- mixOmics::pls(X, Y) # run the method
mixOmics::plotIndiv(pls.result)   # plot the samples
mixOmics::plotVar(pls.result)  
mixOmics::cim(pls.result)


################################################################################
################################################################################
str(vz)
str(fz)

pheatmap::pheatmap(t(vz[mask, ]) %*% fz[mask, ])


################################################################################
################################################################################
################################################################################

################################################################################
# Correlate AA and expression

# Per gene z-Scaled expression
vz <-
  vst.mat |>
  apply(1, scale) |>
  t() |>
  magrittr::set_colnames(colnames(vst.mat))


mask <- intersect(mask, mask2)

crossing(
  AA = colnames(fz),
  lib = colnames(vz)
  # lib = colnames(rz)
) |>
  group_by(AA, lib) |>
  # do(r2 = cor.test(fz[mm, .$AA], rz[mm, .$lib])$estimate) |>
  # do(r2 = cor.test(fz[mask, .$AA], rz[mask, .$lib])$estimate) |>
  do(r2 = cor.test(fz[mask, .$AA], vz[mask, .$lib])$estimate) |>
  # do(r2 = cor.test(freqs.mat[mask, .$AA], vz[mask, .$lib])$estimate) |>
  # do(r2 = cor.test(freqs.mat[mask, .$AA], vz[mask, .$lib])$estimate) |>
  # do(r2 = cor.test(fz[mask, .$AA], vst.mat[mask, .$lib])$estimate) |>
  ungroup() |>
  unnest(r2) |>
  spread(lib, r2) -> dat

dat.mat <- dat |>
  select(- AA) |>
  as.matrix() |>
  magrittr::set_rownames(dat$AA)

################################################################################
# Produce heatmap

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
    set_names(cbPalette[c(1, 6, 2, 7)], .)
)
with(
  meta,
  data.frame('CO2' = as.character(CO2), row.names = sample)
) -> col.df

dat.mat |>
  pheatmap::pheatmap(
    scale = 'none',
    show_colnames = TRUE,
    show_rownames = TRUE,
    display_numbers = TRUE,
    cluster_cols = lib.clust,
    annotation_col = col.df,
    annotation_colors = cl,
    # fontsize = 1,
    number_color = 'black',
    color = colorRampPalette(rev(
      RColorBrewer::brewer.pal(n = 5, name = "RdBu")))(59),
    filename = 'analysis/J_AA-expr-cor.jpeg', width = 8, height = 8
  )
dev.off()
  


################################################################################
sessionInfo()