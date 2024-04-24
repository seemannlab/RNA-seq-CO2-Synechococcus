# Investigate the impact of differential expression on amino acid compositions

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
  # indicate standard deviations
  geom_hline(yintercept = c(-2, 2), linetype = 'dashed') +
  theme(
    legend.position = 'hide',
    axis.text.x = element_text(angle = 60, hjust = 1)
  ) -> p1.z

(p1.total | p1.z) +
  plot_annotation(tag_levels = 'A')

ggsave('analysis/J_freqs-overall.jpeg', width = 14, height = 7, dpi = 400)

################################################################################
# Filter genese if AA does not deviate from average by 2 std, in any AA

mask.aa <-
  fz |>
  abs() |>
  apply( 1, max) %>%
  `>=`(2) |>
#   table()
# FALSE  TRUE 
# 1731  1455 
  which() |>
  names()

################################################################################
# Justify cutoff for lowely expressed, not expression changinge genes

deg |>
  filter(is.de) |>
  group_by(Geneid, name, product, baseMean) |>
  summarize(max.alf = max(abs(log2FoldChange))) |>
  ungroup() -> dat

# top 30% with some rounding on cutoff
dat |>
  select_if(is.numeric) |>
  # map(quantile, probs = seq(0, 1, .1))
  map2(c(.6, .6), ~ quantile(.x, probs = .y)) |>
  map(~ set_names(.x, NULL)) |>
  unlist() |>
  round(1) |>
  as.list() -> xs
xs$baseMean <- round(xs$baseMean, -2) 
# xs

dat |>
  ggplot(aes(baseMean, max.alf)) +
  geom_point(alpha = .5) +
  geom_hline(yintercept = xs$max.alf, color = 'red') +
  geom_vline(xintercept = xs$baseMean, color = 'blue') +
  scale_y_continuous(breaks = 0:8) +
  scale_x_log10(labels = scales::label_comma()) +
  annotate('text', x = 1e4, y = 0, size = 8,
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

# length(mask.aa)
# 1455
# length(mask.deg)
# 440
# length(my.mask)
# 211

################################################################################
# Scatter plot expression to AA freq

# Per gene z-Scaled expression
vz <-
  vst.mat |>
  apply(1, scale) |>
  t() |>
  magrittr::set_colnames(colnames(vst.mat))


vz[my.mask, ] |>
  as_tibble(rownames = 'Geneid') |>
  pivot_longer(- Geneid, names_to = 'lib') |>
  left_join(meta, 'lib') |>
  group_by(Geneid, CO2) |>
  summarize(av = mean(value)) |>
  left_join(
    fz[my.mask, ] |>
      as_tibble(rownames = 'Geneid') |>
      pivot_longer(- Geneid, names_to = 'AA', values_to = 'AA.zfreq'),
    'Geneid',
    relationship = "many-to-many"
  ) |>
  ggplot(aes(av, AA.zfreq)) +
  geom_hex(bins = 20) +
  scale_fill_viridis_c() +
  geom_smooth(method = 'lm', se = FALSE, color = 'red') +
  stat_cor(size = 5, color = 'red') +
  facet_grid(CO2 ~ AA)

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
  sprintf(paste0('%-', sneak.len, 's'), x)
}

dat.mat <- vz[my.mask,]
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
    cluster_cols = heat.cor$tree_col,
    annotation_col = col.df,
    annotation_colors = cl,
    color = colorRampPalette(rev(
      RColorBrewer::brewer.pal(n = 5, name = "RdBu")))(59),
  ) -> heat.expr
dev.off()

################################################################################
# Heatmap AA freq

# apply(abs(dat.cor) >= .1, 1, any) |>
#   which() |>
#   names() -> aas

# dat.mat <- fz[my.mask, aas]
dat.mat <- fz[my.mask, ]
# dat.mat[dat.mat > 4] <- 4
# dat.mat[dat.mat < -4] <- -4


cl <- list(
  'Energy' = c('white', 'black')
)
with(
  aa.cost,
  data.frame('Energy' = cost, row.names = AA)
) -> col.df

# rg <- max(abs(dat.mat)) is 8.6

dat.mat |>
  pheatmap::pheatmap(
    cluster_rows = heat.expr$tree_row,
    cluster_cols = heat.cor$tree_row,
    scale = 'none',
    show_colnames = TRUE,
    show_rownames = FALSE,
    annotation_col = col.df,
    annotation_colors = cl,
    color = colorRampPalette(rev(
      RColorBrewer::brewer.pal(n = 5 , name = "RdBu")))(39),
    # breaks = seq(-rg, rg, length.out = 60),
    # fix to keep extremes darker
    breaks = c(
      -8.7,
      seq(-4, 4, length.out = 38),
      8.7
    )
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
