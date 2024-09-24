# Illustration to better explain the concept of z-scaling to collaborators
# and to explain the "artifact" of averaging per condition

library(tidyverse)
library(ggpubr)
library(patchwork)
library(cowplot)

################################################################################

vst <-
'analysis/D_vst-expression.tsv' |>
read_tsv()

meta <-
  'data/C_meta.tsv' |>
  read_tsv() |>
  mutate_at('CO2', ~ fct_reorder(as.character(.x), as.numeric(.x)))

annot <-
  'data/C_annotation.tsv' |>
  read_tsv()

dat <-
  'analysis/D_normalized-counts.tsv' |>
  read_tsv()

################################################################################
# count matrix

dat.mat <-
  dat |>
  select(- Geneid) |>
  as.matrix() |>
  magrittr::set_rownames(dat$Geneid)


################################################################################
# Heatmap of z-scaling with / without z-scaling

# sample 10 diff genes but from different average expression level
mask <-
  'analysis/K_logFC-vs-30.tsv' |>
  read_tsv() |>
  filter(padj <= 0.001) |>
  arrange(baseMean) |>
  filter(1:n() %% floor(n() / 10) == 0) |>
  pull(Geneid)

x <- log10(dat.mat[mask, ] + 1)

helper <- function(y, path) {
  colnames(x) <-
    x |>
    colnames() |>
    str_remove('Lane.-') |>
    str_remove('_CO2-adaptation')
  rownames(x) <-
    x |>
    rownames() |>
    str_remove('gene-')
  
  pheatmap::pheatmap(
    x,
    scale = y,
    display_numbers = TRUE,
    color = colorRampPalette(rev(
      RColorBrewer::brewer.pal(n = 7, name = "RdBu")))(59),
    filename = path
  )
}

helper('none', 'extra/zscale-noz.jpeg')
helper('row', 'extra/zscale-rowz.jpeg')
dev.off()

# wrap heatmap
p1 <-
  ggdraw() +
  draw_image(
    'extra/zscale-noz.jpeg'
  ) +
  theme_pubr(18) +
  theme(axis.line = element_blank()) +
  ggtitle('log10 expression, no z-scaling')
p2 <-
  ggdraw() +
  draw_image(
    'extra/zscale-rowz.jpeg'
  ) +
  theme_pubr(18) +
  theme(axis.line = element_blank()) +
  ggtitle('per gene (row) z-scaled')

p1 | p2
ggsave('extra/zscle-vis.jpeg', width = 16, height = 9)

################################################################################
