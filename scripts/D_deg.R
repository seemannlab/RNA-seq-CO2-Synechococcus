# Run differential expression analysis for all pairwise condition
# compatison with DESeq2+stageR

library(tidyverse)
library(ggpubr)
library(cowplot)
library(patchwork)

library(DESeq2)
library(stageR)

library(venn)
library(corrplot)

library(conflicted)
conflicts_prefer(dplyr::filter)
conflicts_prefer(dplyr::rename)
conflicts_prefer(purrr::reduce)
conflicts_prefer(dplyr::count)

source('scripts/helper_deg.R')

library(gt)
Sys.setenv(CHROMOTE_CHROME = '/Applications/Brave Browser.app/Contents/MacOS/Brave Browser')

ALPHA <- 0.05

# https://riptutorial.com/r/example/28354/colorblind-friendly-palettes
cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442",
               "#0072B2", "#D55E00", "#CC79A7")

################################################################################
# Load input data

meta <- read_tsv('data/C_meta.tsv') %>%
  mutate_at('CO2', ~ fct_reorder(as.character(.x), as.numeric(.x)))
annot <- read_tsv('data/C_annotation.tsv')
raw.counts <- read_tsv('data/C_raw-counts.tsv')

raw.counts %>%
  select(- Geneid) %>%
  as.matrix() %>%
  magrittr::set_rownames(raw.counts$Geneid) -> raw.counts.mat


################################################################################
# Exclude rRNA

annot %>%
  filter(type != 'rRNA') %>%
  pull(Geneid) -> mask

raw.noribo.mat <- raw.counts.mat[mask, ]

################################################################################
# Run DESeq2+stageR

des <- DESeq2::DESeqDataSetFromMatrix(
  countData = raw.noribo.mat,
  colData = meta,
  design = ~ CO2
)

deg.res <-
  des |>
  helper.pairwise.deg(ALPHA = ALPHA) |>
  helper.stagewise(des = des, ALPHA = ALPHA)

deg.res |>
  write_tsv('analysis/D_stagewise-adjusted-DEGs.tsv')

################################################################################
# Normalized expression levels

des.d <-
  des|>
  DESeq2::DESeq()

des.vst <-
  des.d|>
  DESeq2::vst()

dat.vst <-
  des.vst |>
  SummarizedExperiment::assay()

dat.norm <-
  des.d |>
  DESeq2::counts(normalized = TRUE)

dat.vst |>
  as_tibble(rownames = 'Geneid') |>
  write_tsv('analysis/D_vst-expression.tsv')

dat.norm |>
  as_tibble(rownames = 'Geneid') |>
  write_tsv('analysis/D_normalized-counts.tsv')
  
################################################################################
# PCA


p.pca <- helper.pca(des.vst, sh = 'lane', color_values = cbPalette[c(1, 6, 2, 7)])
p.skree <- helper.skree(des.vst, nudge = 10) + xlim(c(1, 16))


p.pca | p.skree
ggsave('analysis/D_PCA.jpeg', width = 15, height = 8)

################################################################################
# Volcano plot

deg.res |>
  mutate(
    is.de = replace_na(is.de, FALSE),
    is.de = ifelse(is.de, 'yes' , 'no')
  ) |>
  ggplot(aes(log2FoldChange, -log10(padj.stageR), color = is.de)) +
  geom_point(alpha = .7, size = 1) +
  scale_color_manual(name = 'Diff. expressed', values = cbPalette[c(1, 7)]) +
  ylab('-log10 P-Value, stage-wise corrected') +
  facet_wrap(~ test) +
  theme_pubr(18)

ggsave('analysis/D_volcano.jpeg', width = 12, height = 8)

################################################################################
# Heatmap of expression

deg.res %>%
  filter(is.de) %>%
  pull(Geneid) %>%
  unique -> mask

cl <- list(
  'CO2' = meta %>%
    pull(CO2) %>%
    levels %>%
    set_names(cbPalette[c(1, 6, 2, 7)], .)
)

with(
  meta,
  data.frame('CO2' = as.character(CO2), row.names = lib)
) -> cl.df

# Cluster columns ahead of pheatmap to rotate tree nicely
dat <-  dat.vst[mask, ]

# z-scale data ahead of heatmap to keep clustering of rows consistent
dat <-
  dat |>
  apply(1, scale) |>
  t() |>
  magrittr::set_colnames(colnames(dat))

lib.clust <-
  dat |>
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

pheatmap::pheatmap(
  dat,
  scale = 'none',
  cluster_cols = lib.clust,
  show_rownames = FALSE,
  show_colnames = FALSE,
  annotation_col = cl.df,
  annotation_colors = cl,
  color = colorRampPalette(rev(
    RColorBrewer::brewer.pal(n = 7, name = "RdBu")))(59),
  filename = 'analysis/D_heatmap.jpeg'
)
dev.off()




################################################################################
# Correlation of logFCs

deg.res |>
  filter(is.de) |>
  select(Geneid) |>
  unique() |>
  left_join(deg.res, 'Geneid') |>
  select(Geneid, test, log2FoldChange) %>%
  spread(test, log2FoldChange, fill = 0) -> lfc

lfc %>%
  select(- Geneid) %>%
  as.matrix %>%
  magrittr::set_rownames(lfc$Geneid) -> lfc.mat

jpeg('analysis/D_logFC-cor.jpeg', res = 400,
     width = 20, height = 20, units = 'cm')
lfc.mat |>
  cor() |>
  corrplot(
    method = 'square',
    order = 'original',
    type = 'lower',
    diag = FALSE,
    insig='blank',
    addCoef.col ='black',
    number.cex = 0.8,
    col = colorRampPalette(c('#FF6666', '#DDDDDD', '#6666FF'))(50)
  )
dev.off()

################################################################################
###############################################################################
# Overview tables

dat.de <- 
  deg.res |>
  filter(is.de)

dat.total <-
  dat.de |>
  mutate(x = 'DE in any comparison') |>
  select(x, Geneid, type) |>
  unique() |>
  count(x, type) |>
  # sort columns by ammount
  arrange(desc(n))

annot |>
  # highlight Rfam hits with a star
  mutate(is.rfam = str_detect(Geneid, '^Candidate_RF')) |>
  count(Type = type, is.rfam, name = 'Annotated') |>
  arrange(desc(Annotated), Type) |>
  left_join(
    dat.de |>
      select(Geneid, type) |>
      unique() |>
      count(Type = type, name = 'DE'),
    'Type'
  ) |>
  mutate(
    Type = paste0(
      Type,
      ifelse(is.rfam, '*', '')
    ),
    DE = replace_na(DE, 0),
    `%` = DE / Annotated * 100
  ) |>
  select(- is.rfam) -> over.table

over.table |>
  write_tsv('analysis/D_overview.tsv')

###############################################################################
# diff expression by type


dat.de |>
  count(type, test) |>
  mutate_at('type', fct_relevel, dat.total$type) |>
  spread(type, n) |>
  mutate_if(is.numeric, replace_na, 0) %>%
  bind_rows(
    dat.total |>
      rename(test = x) |>
      spread(type, n)
  ) -> type.table


type.table |>
  write_tsv('analysis/D_overview-by-type.tsv')

################################################################################
# No DEGs above logFC

# arrange rows numerically, not by string characters
deg.res |>
  select(test) |>
  unique() |>
  separate(test, c('a', 'b'), sep = '->', remove = FALSE) |>
  mutate_at('b', str_remove, '%.*$') |>
  mutate_at(c('a', 'b'), as.numeric) |>
  arrange(a, b) |>
  drop_na() |>
  pull(test) -> test.ord

list(0, 1, 2, 3, 4, 5) |>
  map(function(i) {
    dat.de |>
      filter(abs(log2FoldChange) >= i) |>
      mutate(x = sprintf('|logFC| ≥ %g', i))
  }) |>
  bind_rows() -> foo

foo |>
  count(x, test) |>
  mutate_at('test', fct_relevel, test.ord) |>
  arrange(test) |>
  spread(x, n, fill = 0) |>
  bind_rows(
    foo |>
      mutate(test = 'DE in any comparison') |>
      select(test, x, Geneid) |>
      unique() |>
      count(test, x) |>
      spread(x, n, fill = 0)
  ) -> lf.table

lf.table |>
  write_tsv('analysis/D_overview-by-logFC.tsv')

################################################################################
################################################################################
# Compile key results into a single pretty overview


# wrap heatmap
p.heat <-
  ggdraw() +
  draw_image(
    'analysis/D_heatmap.jpeg'
  ) +
  theme_pubr(18) +
  theme(axis.line = element_blank())

# prettify tables
lf.table |>
  gt(rowname_col = 'test') |>
  fmt_number(decimals = 0) |>
  data_color(
    rows =  test != 'DE in any comparison',
    palette = "viridis"
  ) |>
  gtsave('tmp-foo.png', zoom = 5)
p.tab <-
  ggdraw() +
  draw_image('tmp-foo.png') +
  theme_pubr(18) +
  theme(axis.line = element_blank())
file.remove('tmp-foo.png')


over.table |>
  gt(rowname_col = 'Type') |>
  fmt_integer() |>
  fmt_percent('%', scale_values = FALSE) |>
  gtsave('tmp-foo.png', zoom = 5)
p.over <-
  ggdraw() +
  draw_image('tmp-foo.png') +
  theme_pubr(18) +
  theme(axis.line = element_blank())
file.remove('tmp-foo.png')

list(
  p.pca + theme(legend.box = 'vertical'),
  wrap_elements(full = p.tab) + theme_pubr(18),
  wrap_elements(full = p.over) + theme_pubr(18),
  wrap_elements(full = p.heat) + theme_pubr(18)
) |>
  # map(~ wrap_elements(full = .x)) |>
  reduce(.f = `+`) +
  plot_layout(
    design = "
    AABBBB
    AADDDD
    CCDDDD
    CCDDDD
    CCDDDD
    "
  ) +
  plot_annotation(tag_levels = 'A')

ggsave('analysis/D_overview.jpeg', width = 11, height = 11, dpi = 600)

################################################################################
################################################################################
sessionInfo()
