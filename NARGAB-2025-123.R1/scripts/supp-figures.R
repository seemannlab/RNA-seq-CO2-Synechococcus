library(tidyverse)
library(ggpubr)
library(patchwork)

library(tidygraph)
library(ggraph)
library(cowplot)

# https://riptutorial.com/r/example/28354/colorblind-friendly-palettes
cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442",
               "#0072B2", "#D55E00", "#CC79A7")

proj.dir <- '/media/dna/projects/rth/co2capture'
my.proj.dir <- file.path(proj.dir, 'continuation-ses/subprojects/RNA-seq_1v15-pct-CO2')

################################################################################
# helpful data for refined figures/tables

meta <-
  file.path(proj.dir, 'subprojects/RNA-seq_1v15-pct-CO2/RNA-seq-CO2-Synechococcus/data/C_meta.tsv') |>
  read_tsv() |>
  mutate_at('CO2', ~ fct_reorder(as.character(.x), as.numeric(.x)))

################################################################################
# some file copies

#file.copy(
#  '~/Local/RNAseq/RNA-seq-CO2-Synechococcus/analysis/D_logFC-cor.jpeg',
#  'files-supplement/fig-logFC-cor.jpeg'
#)

#file.copy(
#  '~/Local/RNAseq/RNA-seq-CO2-Synechococcus/analysis/D_volcano.jpeg',
#  'files-supplement/fig-volcano.jpeg'
#)

#file.copy(
#  '~/Local/RNAseq/RNA-seq-CO2-Synechococcus/analysis/E_treemap-hierarchy.jpeg',
#  'files-supplement/fig-hierarchy.jpeg'
#)


#file.copy(
#  '~/Local/RNAseq/RNA-seq-CO2-Synechococcus/analysis/H_signal-probs.jpeg',
#  'files-supplement/fig-signalP-preds.jpeg'
#)

################################################################################


#file.copy(
#  '~/Local/RNAseq/RNA-seq-CO2-Synechococcus/analysis/I_binomial-model.jpeg',
#  'files-supplement/fig-aa-binomial.jpeg'
#)
#file.copy(
#  '~/Local/RNAseq/RNA-seq-CO2-Synechococcus/analysis/I_correlation-aa-length.jpeg',
#  'files-supplement/fig-aa-cor.jpeg'
#)
#file.copy(
#  '~/Local/RNAseq/RNA-seq-CO2-Synechococcus/analysis/I_negative-binomial-model.jpeg',
#  'files-supplement/fig-aa-nb.jpeg'
#)
#file.copy(
#  '~/Local/RNAseq/RNA-seq-CO2-Synechococcus/analysis/I_volcano-like.jpeg',
#  'files-supplement/fig-aa-volcano.jpeg'
#)
#file.copy(
#  '~/Local/RNAseq/RNA-seq-CO2-Synechococcus/analysis/I_extreme-AA-deg.jpeg',
#  'files-supplement/fig-aa-extreme.jpeg'
#)

################################################################################

#file.copy(
#  '~/Local/RNAseq/RNA-seq-CO2-Synechococcus/analysis/K_logFC-comparisons.jpeg',
#  'files-supplement/fig-focus-cor.jpeg'
#)

#file.copy(
#  '~/Local/RNAseq/RNA-seq-CO2-Synechococcus/analysis/K3_string-cutoff.jpeg',
#  'files-supplement/fig-string-cutoff.jpeg'
#)


file.copy(
  file.path(my.proj.dir, 'RNA-seq-CO2-Synechococcus/analysis/D_tpm-expression-violin-q75.jpeg'),
  file.path(my.proj.dir, 'Manuscript-Overleaf/files-supplement/fig-tpm-expression-violin-q75.jpeg')
)


################################################################################
# Create an illustration of each pair-wise comparison and the screen

xs <-
  meta |>
  pull(CO2) |>
  unique() |>
  sort() |>
  as.character() |>
  paste0('%')

# Regression illustration for screen stage
xs2 <- c(
  'Intercept\n(e.g. 0.04%)',
  xs[- 1]
)

dat <- 
  tibble(
    x = as_factor(xs2),
    y = c(50, 10, 70, 85),
    inter = first(y)
  )
dat |>
  ggplot(aes(x, y)) +
  geom_bar(stat = 'identity', fill = cbPalette[[1]]) +
  geom_hline(
    yintercept = dat |>
      pull(inter) |>
      unique(),
    color = cbPalette[[7]],
    size = 1.5
  ) +
  geom_segment(
    aes(
      x = x,
      y = inter,
      xend = x,
      yend = y,
    ),
    size = 1.5,
    arrow = arrow(angle = 30, length = unit(0.25, "inches"),
                  ends = "both", type = "open"),
    color = cbPalette[[6]],
    # Don't draw for Intercept
    data = slice(dat, -1)
  ) +
  annotate(
    'text',
    x = '4%',
    y = 70,
    label = 'Regression coefficients',
    hjust = .8,
    color = cbPalette[[6]],
    size = 8
  ) +
  xlab(NULL) +
  ylab('Expression level') +
  theme_pubclean(18) -> p1

# Pair-wise

ns <- tibble(name = xs )
es <-
  crossing(from = xs, to = xs) %>%
  mutate(
    x = from |>
      str_remove('%$') |>
      as.numeric(),
    y = to |>
      str_remove('%$') |>
      as.numeric()
  ) |>
  filter(x < y)

grph <- tbl_graph(ns, es)

ggraph(grph, layout = 'linear') +
  geom_edge_arc() +
  geom_node_label(aes(label = name), size = 5) +
  theme_graph(base_size = 18) -> p2
  
p1 + p2 +
  plot_annotation(tag_levels = 'A')

ggsave(
  file.path(proj.dir, 'continuation-ses/subprojects/RNA-seq_1v15-pct-CO2/Manuscript-Overleaf/files-supplement/fig-pairwise.jpeg',
  width = 16, height = 5,
  scale = 1.1,
  dpi = 400
)

################################################################################
################################################################################


################################################################################
# STRING pathways

string.cl2 <- file.path(my.proj.dir, 'RNA-seq-CO2-Synechococcus/string-network/string-30VsOther-cluster2.svg')
svg_image <- image_read_svg(string.cl2, width = 1600)
p.a <- image_ggplot(svg_image)
string.cl4 <- file.path(my.proj.dir, 'RNA-seq-CO2-Synechococcus/string-network/string-30VsOther-cluster4.svg')
svg_image <- image_read_svg(string.cl4, width = 1600)
p.b <- image_ggplot(svg_image)
string.cl7 <- file.path(my.proj.dir, 'RNA-seq-CO2-Synechococcus/string-network/string-30VsOther-cluster7.svg')
svg_image <- image_read_svg(string.cl7, width = 1600)
p.c <- image_ggplot(svg_image)
string.cl6 <- file.path(my.proj.dir, 'RNA-seq-CO2-Synechococcus/string-network/string-30VsOther-cluster6.svg')
svg_image <- image_read_svg(string.cl6, width = 1600)
p.d <- image_ggplot(svg_image)
string.cl9 <- file.path(my.proj.dir, 'RNA-seq-CO2-Synechococcus/string-network/string-30VsOther-cluster9.svg')
svg_image <- image_read_svg(string.cl9, width = 1600)
p.e <- image_ggplot(svg_image)
string.cl14 <- file.path(my.proj.dir, 'RNA-seq-CO2-Synechococcus/string-network/string-30VsOther-cluster14.svg')
svg_image <- image_read_svg(string.cl14, width = 1600)
p.f <- image_ggplot(svg_image)
string.cl15 <- file.path(my.proj.dir, 'RNA-seq-CO2-Synechococcus/string-network/string-30VsOther-cluster15.svg')
svg_image <- image_read_svg(string.cl15, width = 1600)
p.g <- image_ggplot(svg_image)
string.cl8 <- file.path(my.proj.dir, 'RNA-seq-CO2-Synechococcus/string-network/string-30VsOther-cluster8.svg')
svg_image <- image_read_svg(string.cl8, width = 1600)
p.h <- image_ggplot(svg_image)

design <- "AAAAABBBBCCC
           AAAAABBBBCCC
           AAAAABBBBCCC
           AAAAABBBBCCC
           AAAAABBBBCCC
           DDDEEFFGGHHH
           DDDEEFFGGHHH
           DDDEEFFGGHHH"
p.a + p.b + p.c + p.d + p.e + p.f + p.g + p.h +
  plot_layout(design = design) + plot_annotation(tag_levels = 'A') & theme(plot.tag = element_text(size = 40))
ggsave(file.path(proj.dir, 'continuation-ses/subprojects/RNA-seq_1v15-pct-CO2/Manuscript-Overleaf/files-supplement/fig-string-more-clusters.pdf'), width = 20, height = 12)


################################################################################
# jbrowser images of CRSs, Rfam, and orphan RNAs

# colors used in JBrowser
my.cl <- list(
  '0.04' = c(160, 160, 160),
  '4' = c(14, 121, 182),
  '8' = c(231, 162, 9),
  '30' = c(214, 99, 10)
) |>
  map(~ invoke(rgb, .x, maxColorValue = 255))

# custom CO2 condition color legend
p.leg <-
  tibble(
    co2 = names(my.cl) |> fct_inorder(),
  ) |>
  ggplot(aes(co2, 1, fill = co2)) +
  scale_fill_manual(values = my.cl, name = bquote('%'~CO[2])) +
  geom_tile() +
  guides(fill = guide_legend(nrow = 1)) +
  theme_pubr(18)

p.leg <- ggpubr::get_legend(p.leg)

png.dir <- file.path(my.proj.dir, "RNA-seq-CO2-Synechococcus/fig-orphan-crs")

#annotated sRNAs (Rfam)
rna.1 <- file.path(png.dir, 'jbrowse-tmrna.png')
p.rna.1 <- ggdraw() + draw_image(magick::image_read(rna.1, density = 100)) + theme(plot.margin = margin(t=.5, l=.5, r=.5, b=.5, unit = "cm"))
rna.2 <- file.path(png.dir, 'jbrowse-RNaseP.png')
p.rna.2 <- ggdraw() + draw_image(magick::image_read(rna.2, density = 100)) + theme(plot.margin = margin(t=.5, l=.5, r=.5, b=.5, unit = "cm"))
rna.3 <- file.path(png.dir, 'jbrowse-6Srna.png')
p.rna.3 <- ggdraw() + draw_image(magick::image_read(rna.3, density = 100)) + theme(plot.margin = margin(t=.5, l=.5, r=.5, b=.5, unit = "cm"))
rna.4 <- file.path(png.dir, 'jbrowse-SRP.png')
p.rna.4 <- ggdraw() + draw_image(magick::image_read(rna.4, density = 100)) + theme(plot.margin = margin(t=.5, l=.5, r=.5, b=.5, unit = "cm"))
rna.5 <- file.path(png.dir, 'jbrowse-PsrR1.png')
p.rna.5 <- ggdraw() + draw_image(magick::image_read(rna.5, density = 100)) + theme(plot.margin = margin(t=.5, l=.5, r=.5, b=.5, unit = "cm"))

#annotated riboswitches
rna.6 <- file.path(png.dir, 'jbrowse-TPP-riboswitch.png')
p.rna.6 <- ggdraw() + draw_image(magick::image_read(rna.6, density = 100)) + theme(plot.margin = margin(t=.5, l=.5, r=.5, b=.5, unit = "cm"))
rna.7 <- file.path(png.dir, 'jbrowse-ydaO-yuaA-riboswitch.png')
p.rna.7 <- ggdraw() + draw_image(magick::image_read(rna.7, density = 100)) + theme(plot.margin = margin(t=.5, l=.5, r=.5, b=.5, unit = "cm"))

#CRSs
rna.8.1 <- file.path(png.dir, 'jbrowse-K02110.png')
p.rna.8.1 <- ggdraw() + draw_image(magick::image_read(rna.8, density = 100)) + theme(plot.margin = margin(t=.5, l=.5, r=.5, b=.5, unit = "cm"))
rna.8.2 <- file.path(png.dir, 'K02110_downstream.h1_5.png')
p.rna.8.2 <- ggdraw() + draw_image(magick::image_read(rna.8.2, density = 100)) + theme(plot.margin = margin(t=.5, l=.5, r=.5, b=.5, unit = "cm"))
p.rna.8 <- plot_grid(p.rna.8.1,p.rna.8.2, ncol=2, rel_widths = c(10,1))
rna.9.1 <- file.path(png.dir, 'jbrowse-K02094.png')
p.rna.9.1 <- ggdraw() + draw_image(magick::image_read(rna.9, density = 100)) + theme(plot.margin = margin(t=.5, l=.5, r=.5, b=.5, unit = "cm"))
rna.9.2 <- file.path(png.dir, 'K02094_upstream.h1_2.h1_4.png')
p.rna.9.2 <- ggdraw() + draw_image(magick::image_read(rna.9.2, density = 100)) + theme(plot.margin = margin(t=.5, l=.5, r=.5, b=.5, unit = "cm"))
p.rna.9 <- plot_grid(p.rna.9.1,p.rna.9.2, ncol=2, rel_widths = c(3,1.2))
rna.10.1 <- file.path(png.dir, 'jbrowse-K07086.png')
p.rna.10.1 <- ggdraw() + draw_image(magick::image_read(rna.10, density = 100)) + theme(plot.margin = margin(t=.5, l=.5, r=.5, b=.5, unit = "cm"))
rna.10.2 <- file.path(png.dir,'K07086_upstream.h2_3.png')
p.rna.10.2 <- ggdraw() + draw_image(magick::image_read(rna.10.2, density = 100)) + theme(plot.margin = margin(t=.5, l=.5, r=.5, b=.5, unit = "cm"))
p.rna.10 <- plot_grid(p.rna.10.1,p.rna.10.2, ncol=2, rel_widths = c(3,1))
#rna.11 <- file.path(png.dir, 'jbrowse-K21903.png')
#p.rna.11 <- ggdraw() + draw_image(magick::image_read(rna.11, density = 100)) + theme(plot.margin = margin(t=.5, l=.5, r=.5, b=.5, unit = "cm"))

#orphan sRNAs
npc01 <- file.path(png.dir, 'jbrowse-RS02500_RS02510.png')
p.npc01 <- ggdraw() + draw_image(magick::image_read(npc01, density = 100)) + theme(plot.margin = margin(t=.5, l=.5, r=.5, b=.5, unit = "cm"))
npc02 <- file.path(png.dir, 'jbrowse-RS02980_RS02985.png')
p.npc02 <- ggdraw() + draw_image(magick::image_read(npc02, density = 100)) + theme(plot.margin = margin(t=.5, l=.5, r=.5, b=.5, unit = "cm"))
npc03 <- file.path(png.dir, 'jbrowse-RS06725_RS06730.png')
p.npc03 <- ggdraw() + draw_image(magick::image_read(npc03, density = 100)) + theme(plot.margin = margin(t=.5, l=.5, r=.5, b=.5, unit = "cm"))
npc04 <- file.path(png.dir, 'jbrowse-RS09150_RS09155.png')
p.npc04 <- ggdraw() + draw_image(magick::image_read(npc04, density = 100)) + theme(plot.margin = margin(t=.5, l=.5, r=.5, b=.5, unit = "cm"))
npc05 <- file.path(png.dir, 'jbrowse-RS10475_RS10480.png')
p.npc05 <- ggdraw() + draw_image(magick::image_read(npc05, density = 100)) + theme(plot.margin = margin(t=.5, l=.5, r=.5, b=.5, unit = "cm"))
npc06 <- file.path(png.dir, 'jbrowse-RS11335_RS11340.png')
p.npc06 <- ggdraw() + draw_image(magick::image_read(npc06, density = 100)) + theme(plot.margin = margin(t=.5, l=.5, r=.5, b=.5, unit = "cm"))
npc07 <- file.path(png.dir, 'jbrowse-RS13490_RS13495.png')
p.npc07 <- ggdraw() + draw_image(magick::image_read(npc07, density = 100)) + theme(plot.margin = margin(t=.5, l=.5, r=.5, b=.5, unit = "cm"))
npc08 <- file.path(png.dir, 'jbrowse-RS13890_RS13895.png')
p.npc08 <- ggdraw() + draw_image(magick::image_read(npc08, density = 100)) + theme(plot.margin = margin(t=.5, l=.5, r=.5, b=.5, unit = "cm"))
npc09 <- file.path(png.dir, 'jbrowse-RS14650_RS14655.png')
p.npc09 <- ggdraw() + draw_image(magick::image_read(npc09, density = 100)) + theme(plot.margin = margin(t=.5, l=.5, r=.5, b=.5, unit = "cm"))

#5'-ext+up+switch
rna.21 <- file.path(png.dir, 'jbrowse-RS03930_RS03935.png')
p.rna.21 <- ggdraw() + draw_image(magick::image_read(rna.21, density = 100)) + theme(plot.margin = margin(t=.5, l=.5, r=.5, b=.5, unit = "cm"))
rna.22 <- file.path(png.dir, 'jbrowse-RS01315_RS01320.png')
p.rna.22 <- ggdraw() + draw_image(magick::image_read(rna.22, density = 100)) + theme(plot.margin = margin(t=.5, l=.5, r=.5, b=.5, unit = "cm"))

#putative orphan peptides
rna.23 <- file.path(png.dir, 'jbrowse-RS03600_RS03605.png')
p.rna.23 <- ggdraw() + draw_image(magick::image_read(rna.23, density = 100)) + theme(plot.margin = margin(t=.5, l=.5, r=.5, b=.5, unit = "cm"))

plot_grid(
  p.leg,
  ggarrange(p.rna.1, p.rna.2, p.rna.3, p.rna.4, p.rna.5,
          labels = c("A", "B", "C", "D", "E"),
          ncol = 2, nrow = 3),
  ncol = 1,
  rel_heights = c(1, 10)
)
ggsave(
  file.path(proj.dir, 'continuation-ses/subprojects/RNA-seq_1v15-pct-CO2/Manuscript-Overleaf/files-supplement/fig-rfam-srna.jpeg'),
      width = 8, height = 7,
      scale = 1.1,
      dpi = 400
)

plot_grid(
  p.leg,
  ggarrange(p.rna.6, p.rna.7,
            labels = c("A", "B"),
            ncol = 1, nrow = 2),
  ncol = 1,
  rel_heights = c(1, 10)
)
ggsave(
  file.path(proj.dir, 'continuation-ses/subprojects/RNA-seq_1v15-pct-CO2/Manuscript-Overleaf/files-supplement/fig-riboswitch.jpeg'),
  width = 8, height = 7,
  scale = 1.1,
  dpi = 400
)

plot_grid(
  p.leg,
  ggarrange(p.rna.8, p.rna.9, p.rna.10,
            labels = c("A", "B", "C"),
            ncol = 1, nrow = 3),
  ncol = 1,
  rel_heights = c(1, 15)
)
ggsave(
  file.path(proj.dir, 'continuation-ses/subprojects/RNA-seq_1v15-pct-CO2/Manuscript-Overleaf/files-supplement/fig-crs.jpeg'),
  width = 8, height = 11,
  scale = 1.1,
  dpi = 400
)

plot_grid(
  p.leg,
  ggarrange(p.npc01, p.npc02, p.npc03, p.npc04, p.npc05,p.npc06, p.npc07, p.npc08, p.npc09,
            labels = c("A", "B", "C", "D", "E", "F", "G", "H", "I"),
            ncol = 2, nrow = 5),
  ncol = 1,
  rel_heights = c(1, 18)
)
ggsave(
  file.path(proj.dir, 'continuation-ses/subprojects/RNA-seq_1v15-pct-CO2/Manuscript-Overleaf/files-supplement/fig-orphan-rna.jpeg'),
  width = 8, height = 11,
  scale = 1.1,
  dpi = 400
)

plot_grid(
  p.leg,
  ggarrange(p.rna.21, p.rna.22,
            labels = c("A", "B"),
            ncol = 1, nrow = 2),
  ncol = 1,
  rel_heights = c(1, 10)
)
ggsave(
  file.path(proj.dir, 'continuation-ses/subprojects/RNA-seq_1v15-pct-CO2/Manuscript-Overleaf/files-supplement/fig-5ext-up-switch.jpeg'),
  width = 8, height = 7,
  scale = 1.1,
  dpi = 400
)

plot_grid(
  p.leg,
  p.rna.23,
  ncol = 1,
  rel_heights = c(1, 10)
)
ggsave(
  file.path(proj.dir, 'continuation-ses/subprojects/RNA-seq_1v15-pct-CO2/Manuscript-Overleaf/files-supplement/fig-orphan-peptide.jpeg'),
  width = 8, height = 3.8,
  scale = 1.2,
  dpi = 400
)

