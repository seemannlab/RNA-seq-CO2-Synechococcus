# Add expression levels for the genes in the pretty carotenoid pathway

library(tidyverse)
library(cowplot)

################################################################################

meta <-
  'data/C_meta.tsv' |>
  read_tsv() |>
  mutate(CO2 = CO2 |>
           as_factor() |>
           fct_reorder(CO2))

vst <-
  'analysis/D_vst-expression.tsv' |>
  read_tsv()

annot <-
  'data/C_annotation.tsv' |>
  read_tsv()

################################################################################

pathway <-
  ggdraw() +
  draw_image(magick::image_read(
    'fig-carotenoid-pathway/carotenoid.jpg',
    density = 500
  ))


################################################################################
# coordiantes and genes as in idraw/svg from top left cornder in pixels

genes <- tribble(
  ~ gene, ~ locus, ~x, ~y,
  'crtE', 'A1085', 1220, 170,
  'crtB', 'A1936', 1220, 310,
  'crtP', 'A1935', 1220, 450,
  'crtQ', 'A0529', 1220, 590,
  'crtH', 'A1829', 1220, 730,
  
  'cruF', 'A2032', 610, 520,
  
  'cruA', 'A2153', 600, 1090,
  'cruP', 'A0043', 600, 1140,
  
  'crtR', 'A0915', 600, 1320,
  
  'cruG', 'A2031', 30, 1440,
  
  'cruE', 'A1248', 1320, 920,
  'crtW', 'A2809', 1220, 1430,
  'cruH', 'A2246', 2000, 1400,
) |>
  # convert to coordinates in ggplot
  # svg is
  # height: 1713 px
  # width: 2271 px
  mutate(
    x = x / 2271,
    y = 1 - y / 1713,
    locus = paste0('SYNPCC7002_', locus)
  )

# overall positions seem to fit
# pathway +
#   geom_text(aes(x, y, label = gene), data = genes)

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
  left_join(meta, 'lib') |>
  group_by(Geneid, CO2) |>
  summarize(avg = mean(z.expression)) |>
  ungroup() |>
  left_join(annot, 'Geneid') |>
  inner_join(genes, c('old_locus_tag' = 'locus'))

################################################################################
################################################################################

# create a small heatmap for a single gene

small.heatmap <- function(dat, width = 0.1, height = 0.035) {
  # dat <- filter(z.expr, Geneid == 'gene-SYNPCC7002_RS00220')
  p <-
    dat |>
    ggplot(aes(CO2, 1, fill = avg)) +
    geom_tile(color = 'black', linewidth = 3) +
    coord_fixed() +
    scale_fill_gradientn(colors = RColorBrewer::brewer.pal(5, 'RdBu') |> rev(),
                         # make sure all have the same scale
                         limits = c(-2, 2)) +
    theme_void() +
    theme(legend.position = 'hide')
  # helper annotation to add to main plot
  dat |>
    slice(1) |>
    with(annotation_custom(
      ggplotGrob(p), 
      xmin = x,   xmax = x + width,
      ymin = y - height, ymax = y
    ))
}
  
heats <-
  z.expr |>
  group_by(old_locus_tag) |>
  do(i = small.heatmap(.)) |>
  pull(i)
  
  
# create legend
p.legend <-
  z.expr |>
  # arbitrarily choose first one for color bar
  filter(Geneid == first(Geneid)) |>
  ggplot(aes(CO2, 1, fill = avg, label = CO2)) +
  geom_tile() +
  coord_fixed() +
  scale_fill_gradientn(colors = RColorBrewer::brewer.pal(5, 'RdBu') |> rev(),
                       name = 'Average z-scaled expression',
                       limits = c(-2, 2)) +
  guides(fill = guide_colorbar(title.position = 'top', barwidth = unit(25, 'cm'))) +
  ggpubr::theme_pubclean(62)  +
  theme(
    axis.title.y = element_blank(),
    axis.text.y =  element_blank(),
    axis.ticks = element_blank(),
    panel.grid = element_blank(),
    panel.grid.major.y = element_blank()
  ) +
  geom_tile(fill = I('white'), color = 'black', linewidth = 3) +
  geom_text(size = 18) +
  theme(axis.text.x = element_blank())
  
# combine all plots into one
p.all <-
  heats |>
  reduce(`+`, .init = pathway) +
  annotation_custom(
    ggplotGrob(p.legend),
    xmin = 0.7, xmax = 0.95,
    ymin = 0.7, ymax = 0.9
  )

ggsave(
  'fig-carotenoid-pathway/carotenoid-expression.jpg', plot = p.all,
  height = 11896,
  width = 15771,
  limitsize = FALSE,
  units = 'px'
)
  
