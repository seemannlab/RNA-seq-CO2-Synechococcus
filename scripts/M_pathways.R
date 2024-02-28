#!/usr/bin/env Rscript

# Plot KEGG pathways with added expression heatmaps
library(tidyverse)
library(ggkegg)

# Pathways of interest
# Glyoxylate and dicarboxylate metabolism - Picosynechococcus sp. PCC 7002 
# https://www.kegg.jp/pathway/syp00630
# Carotenoid biosynthesis 
# https://www.kegg.jp/pathway/syp00906+SYNPCC7002_A1248
# Biotin metabolism
# https://www.kegg.jp/pathway/syp00780

################################################################################

annot <-
  'data/C_annotation.tsv' |>
  read_tsv()

vst <-
  'analysis/E_vst.tsv' |>
  read_tsv()

meta <-
  'data/C_meta.tsv' |>
  read_tsv() |>
  mutate_at('condition', ~ fct_reorder(as.character(.x), as.numeric(.x)))

################################################################################
# z-scale expression

vst.mat <-
  vst |>
  select(- Geneid) |>
  as.matrix() |>
  magrittr::set_rownames(vst$Geneid)

z.mat <-
  vst.mat |>
  apply(1, scale) |>
  t() |> 
  magrittr::set_colnames(colnames(vst.mat))

z.expr <-
  z.mat |>
  as_tibble(rownames = 'Geneid') |>
  pivot_longer(- Geneid, names_to = 'lib', values_to = 'z.expression') |>
  left_join(meta, 'lib')

################################################################################
my.path <- function(x, legend.pos) {
  # Memo: Legend.pos is list with xmin/xmax/ymin/ymax
  # x <- 'syp00630'
  # Load pathway
  path <-
    x |>
    pathway()
  
  # extract genes from pathway
  path.genes <-
    path |>
    activate(nodes) |>
    as_tibble() |>
    filter(type == 'gene')
  # join to expression
  path.expr <-
    path.genes |>
    select(name) |>
    unique() |>
    # match locus to Geneid
    mutate(locus = str_remove_all(name, 'syp:')) |>
    separate_rows(locus, sep = ' ') |>
    left_join(annot, c('locus' = 'old_locus_tag')) |>
    select(name = name.x, Geneid) |>
    # add z-scaled expression
    left_join(z.expr, 'Geneid', relationship = "many-to-many") |>
    # average per 'enzyme' NOT gene (could be multiple)
    group_by(name, condition) |>
    summarize(avg = mean(z.expression, na.rm = TRUE)) |>
    ungroup() |>
    # left_join(path.genes, 'name', relationship = "many-to-many") |>
    # ensure all values fit the ±2 scale
    mutate(avg = pmin(pmax(avg, -2), 2))
  
  # base KEGG path viz
  p.path <-
    path |>
    ggraph(layout="manual", x=x, y=y) +
    # highlight genes in organism
    ggfx::with_outer_glow(
      geom_node_rect(aes(filter = (type == 'gene')), fill = 'yellow'),
      colour="yellow",expand=5
    ) +
    overlay_raw_map() +
    theme_void()
  # p.path
  
  # add expression for genes into the base visualization
  offset <- c(2, - 13)
  helper <- function(gene_row) {
    p <-
      gene_row |>
      as_tibble() |>
      select(name) |>
      left_join(path.expr, 'name', relationship = 'one-to-many') |>
      ggplot(aes(condition, 1, fill = avg)) +
      geom_tile(color = 'black') +
      coord_fixed() +
      scale_fill_gradientn(colors = RColorBrewer::brewer.pal(5, 'RdBu') |> rev(),
                           # make sure all have the same scale
                           limits = c(-2, 2)) +
      theme_void() +
      theme(legend.position = 'hide')
    # helper annotation to add to main plot
    gene_row |>
      with(annotation_custom(
        ggplotGrob(p), 
        xmin = xmin - offset[1], xmax = xmax + offset[1],
        ymin = ymin + offset[2], ymax = ymax + offset[2]
      ))
  }
  
  gene.plots <-
    path.genes |>
    rowwise() |>
    do(i = helper(.)) |>
    pull(i)
  
  
  # create legend
  p.legend <-
    path.expr |>
    # arbitrarily choose first one for color bar
    filter(name == first(name)) |>
    ggplot(aes(condition, 1, fill = avg, label = condition)) +
    geom_tile() +
    coord_fixed() +
    scale_fill_gradientn(colors = RColorBrewer::brewer.pal(5, 'RdBu') |> rev(),
                         name = 'Average z-scaled expression',
                         limits = c(-2, 2)) +
    guides(fill = guide_colorbar(title.position = 'top', barwidth = unit(5, 'cm'))) +
    ggpubr::theme_pubclean(8)  +
    theme(
      axis.title.y = element_blank(),
      axis.text.y =  element_blank(),
      axis.ticks = element_blank(),
      panel.grid = element_blank(),
      panel.grid.major.y = element_blank()
    ) +
    geom_tile(fill = I('white'), color = 'black') +
    geom_text(size = 3) +
    theme(axis.text.x = element_blank())
  
  # combine all plots into one
  gene.plots |>
    reduce(`+`, .init = p.path) +
    with(
      legend.pos,
      annotation_custom(
        ggplotGrob(p.legend),
        xmin = xmin, xmax = xmax,
        ymin = ymin, ymax = ymax
      )
    )
}
  
################################################################################
  
# Glyoxylate and dicarboxylate metabolism
foo <- my.path(
  'syp00630',
  legend.pos = list(
    xmin = 1000, xmax = 1300,
    ymin = -850, ymax = -1050
  )
)
# remove ugly thick border
foo + theme(plot.margin = unit(c(-1, -2, -1, -2), "cm"))

ggsave('analysis/M_photorespiration.jpeg',
       width = 12.5, height = 10, dpi = 400)

########################################

# Carotenoid biosynthesis 
bar <- my.path(
  'syp00906',
  legend.pos = list(
    xmin = 260, xmax = 560,
    ymin = -5, ymax = -160
  )
)
# remove ugly thick border
bar + theme(plot.margin = unit(c(-1, -2, -1, -2), "cm"))
ggsave('analysis/M_carotenoid.jpeg',
       width = 14, height = 10, dpi = 400)

########################################

# Biotin metabolism
baz <-  my.path(
  'syp00780',
  legend.pos = list(
    xmin = 500, xmax = 700,
    ymin = 0, ymax = -200
  )
)
  
# remove ugly thick border
baz + theme(plot.margin = unit(c(-1, -2, -1, -2), "cm"))
ggsave('analysis/M_biotin.jpeg',
       width = 8, height = 10, dpi = 400)
