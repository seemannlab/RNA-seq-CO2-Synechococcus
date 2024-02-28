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
  read_tsv()

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
  theme(
    legend.position = 'hide',
    axis.text.x = element_text(angle = 60, hjust = 1)
  ) -> p1.z

(p1.total | p1.z) +
  plot_annotation(tag_levels = 'A')

ggsave('analysis/J_freqs-overall.jpeg', width = 14, height = 7, dpi = 400)

################################################################################
# Correlate AA and expression

# Per gene z-Scaled expression
vz <-
  vst.mat |>
  apply(1, scale) |>
  t() |>
  magrittr::set_colnames(colnames(vst.mat))


# Focus on DEGs
deg |>
  # filter(is.de, abs(log2FoldChange) >= 1) |>
  filter(is.de) |>
  pull(Geneid) |>
  unique() |>
  intersect(rownames(fz)) -> mask

crossing(
  AA = colnames(fz),
  lib = colnames(vz)
) |>
  group_by(AA, lib) |>
  do(r2 = cor.test(fz[mask, .$AA], vz[mask, .$lib])$estimate) |>
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