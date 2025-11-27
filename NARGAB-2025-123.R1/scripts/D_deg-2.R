# Run differential expression analysis for all pairwise conditions and 30% versus rest

# Purpose: Expression analysis, but focused on differences to 30% using DESeq2 with individual contrast averages

library(tidyverse)
library(gt)

library(cowplot)
library(patchwork)
library(ggpubr)

library(DESeq2)

library(venn)
library(corrplot)

library(enrichplot)
library(clusterProfiler)

library(conflicted)
conflicts_prefer(dplyr::filter)
conflicts_prefer(dplyr::rename)
conflicts_prefer(purrr::reduce)
conflicts_prefer(dplyr::count)

input.dir <- "/media/dna/projects/rth/co2capture/subprojects/RNA-seq_1v15-pct-CO2/RNA-seq-CO2-Synechococcus/"
output.dir <- "/media/dna/projects/rth/co2capture/continuation-ses/subprojects/RNA-seq_1v15-pct-CO2/RNA-seq-CO2-Synechococcus/"
source(paste0(output.dir, 'scripts/helper_deg.R'))

# https://riptutorial.com/r/example/28354/colorblind-friendly-palettes
cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442",
               "#0072B2", "#D55E00", "#CC79A7")

################################################################################
# Load input data

annot <-
  paste0(input.dir, 'data/C_annotation.tsv') |>
  read_tsv()

meta <-
  paste0(input.dir, 'data/C_meta.tsv') |>
  read_tsv() |>
  mutate_at('CO2', ~ fct_reorder(as.character(.x), as.numeric(.x)))

raw.counts <-
  paste0(input.dir, 'data/C_raw-counts.tsv') |>
  read_tsv()


################################################################################
# Build matrices

raw.counts.mat <-
  raw.counts %>%
  select(- Geneid) %>%
  as.matrix() %>%
  magrittr::set_rownames(raw.counts$Geneid)


################################################################################
# Exclude rRNA, CRSs and CRISPR-DR33 (RF01343)

annot %>%
  filter(type != 'rRNA' & type != 'CRS' & type!= 'CRISPR-DR33') %>%
  pull(Geneid) -> mask

raw.noribo.mat <- raw.counts.mat[mask, ]


################################################################################
# Transcripts Per Kilobase Million (TPM) and highly expressed genes

annot <- annot %>% mutate(length = end - start + 1)
annot.l <- annot %>% column_to_rownames("Geneid")
annot.l <- annot.l[mask,] %>% select(length)

raw.noribo.mat |>
  as_tibble(rownames = 'Geneid') |>
  left_join(annot, 'Geneid') |>
  select(
    Geneid,
    'S1' = `CO2-0-04percent_Lane1-S1_CO2-adaptation`,
    'S2' = `CO2-0-04percent_Lane2-S2_CO2-adaptation`,
    'S3' = `CO2-0-04percent_Lane2-S3_CO2-adaptation`,
    'S4' = `CO2-0-04percent_Lane2-S4_CO2-adaptation`,
    'S5' = `CO2-4percent_Lane2-S5_CO2-adaptation`,
    'S6' = `CO2-4percent_Lane1-S6_CO2-adaptation`,
    'S7' = `CO2-4percent_Lane1-S7_CO2-adaptation`,
    'S8' = `CO2-4percent_Lane2-S8_CO2-adaptation`,
    'S9' = `CO2-8percent_Lane2-S9_CO2-adaptation`,
    'S10' = `CO2-8percent_Lane2-S10_CO2-adaptation`,
    'S11' = `CO2-8percent_Lane2-S11_CO2-adaptation`,
    'S12' = `CO2-8percent_Lane1-S12_CO2-adaptation`,
    'S13' = `CO2-30percent_Lane2-S13_CO2-adaptation`,
    'S14' = `CO2-30percent_Lane1-S14_CO2-adaptation`,
    'S15' = `CO2-30percent_Lane2-S15_CO2-adaptation`,
    'S16' = `CO2-30percent_Lane1-S16_CO2-adaptation`,
    length
  ) |>
  write_tsv(paste0(output.dir, 'analysis/D_raw-expression-noribo.tsv'))

library(DGEobj.utils)
#see https://haroldpimentel.wordpress.com/2014/05/08/what-the-fpkm-a-review-rna-seq-expression-units/
#countToTpm <- function(counts, effLen)
#{
#  rate <- log(counts) - log(effLen)
#  denom <- log(sum(exp(rate)))
#  exp(rate - denom + log(1e6))
#}
tpm <- convertCounts(as.matrix(raw.noribo.mat), unit = "TPM", geneLength = annot.l[,1])
tpm <- tpm %>% as.data.frame() %>%
  rownames_to_column("rowname") %>%
  rowwise() %>% mutate(min_tpm = min(c_across(-rowname)), max_tpm = max(c_across(-rowname))) %>%
  column_to_rownames("rowname")

tbl.tpm <-
  tpm |>
  dplyr::rename(
    'S1' = `CO2-0-04percent_Lane1-S1_CO2-adaptation`,
    'S2' = `CO2-0-04percent_Lane2-S2_CO2-adaptation`,
    'S3' = `CO2-0-04percent_Lane2-S3_CO2-adaptation`,
    'S4' = `CO2-0-04percent_Lane2-S4_CO2-adaptation`,
    'S14' = `CO2-30percent_Lane1-S14_CO2-adaptation`,
    'S16' = `CO2-30percent_Lane1-S16_CO2-adaptation`,
    'S13' = `CO2-30percent_Lane2-S13_CO2-adaptation`,
    'S15' = `CO2-30percent_Lane2-S15_CO2-adaptation`,
    'S6' = `CO2-4percent_Lane1-S6_CO2-adaptation`,
    'S7' = `CO2-4percent_Lane1-S7_CO2-adaptation`,
    'S5' = `CO2-4percent_Lane2-S5_CO2-adaptation`,
    'S8' = `CO2-4percent_Lane2-S8_CO2-adaptation`,
    'S12' = `CO2-8percent_Lane1-S12_CO2-adaptation`,
    'S10' = `CO2-8percent_Lane2-S10_CO2-adaptation`,
    'S11' = `CO2-8percent_Lane2-S11_CO2-adaptation`,
    'S9' = `CO2-8percent_Lane2-S9_CO2-adaptation`,
    'min TPM' = 'min_tpm',
    'max TPM' = 'max_tpm'
  ) |>
  relocate(S5, .before = S6) |> relocate(S9, .before = S10) |> relocate(S12, .after = S11) |> relocate(S13, .after = S12) |> relocate(S14, .after = S13) |> relocate(S15, .after = S14) |> relocate(S16, .after = S15)

#get statistics to highly expressed genes
#gather data into long format
tpm.long <- tbl.tpm %>% select(1:(ncol(tbl.tpm)-2)) %>% pivot_longer(cols = everything(), names_to = "Sample", values_to = "TPM")
sample.order <- c("S1", "S2", "S3", "S4", "S5", "S6", "S7", "S8", "S9", "S10", "S11", "S12", "S13", "S14", "S15", "S16")
tpm.long$Sample <- factor(tpm.long$Sample, levels = sample.order)

#violine plot with marked quantiles

ggplot(tpm.long, aes(x = Sample, y = TPM)) + geom_violin(draw_quantiles = c(.75)) + labs(x = "Samples", y = "TPM (log10 scaled)") + scale_y_log10() + theme_pubr(18)

ggsave(paste0(output.dir, 'analysis/D_tpm-expression-violin-q75.jpeg'), width = 12, height = 8)

#calculate quantiles for each column
quantiles <- tpm.long %>% group_by(Sample) %>% summarize(q25 = quantile(TPM, .25), q50 = quantile(TPM, .5), q75 = quantile(TPM, .75), q90 = quantile(TPM, .9))
thres.high.expr <- min(quantiles$q75)
tbl.tpm <- tbl.tpm %>% mutate(`high expr` = `max TPM` > min(quantiles$q75))
tbl.tpm %>% group_by(`high expr`) %>% summarise(n = tidygraph::n())

tbl.tpm |>
  as_tibble(rownames = 'Geneid') |>
  write_tsv(paste0(output.dir, 'analysis/D_tpm-expression.tsv'))


################################################################################
# Run DESeq2 for all pairwise comparisons to check how much 4% and 8% CO2 differ

des <- DESeq2::DESeqDataSetFromMatrix(
  countData = raw.noribo.mat,
  colData = meta,
  design = ~ CO2 + 0
)
des <- DESeq2::DESeq(des)

deg.res.pw <-
  des |>
  helper.pairwise.deg(ALPHA = 0.05) |>
  rename("Geneid" = row) |>
  mutate("is.de" = ifelse(padj <= 0.05, TRUE, FALSE))

deg.res.pw %>% filter(test == "4->8% CO2" & is.de & abs(log2FoldChange)>=1) %>% nrow() #-> 1 DEG


################################################################################
# Run DESeq2 for all pairwise comparisons and 30% VS Rest after combining optimal growth conditions (4 and 8%)

meta2 <-
  meta |>
  mutate(CO2 = case_when(CO2 %in% c(4,8) ~ '4.8',
                         CO2 == 0.04 ~ '0.04',
                         CO2 == 30 ~ '30') |>
           fct_relevel(c('0.04', '4.8', '30')))

des <- DESeq2::DESeqDataSetFromMatrix(
  countData = raw.noribo.mat,
  colData = meta2,
  design = ~ CO2 + 0
)
des <- DESeq2::DESeq(des)

#all pairwise comparisons
deg.res <-
  des |>
  helper.pairwise.deg(ALPHA = 0.05) |>
  rename(Geneid = row) |>
  mutate("is.de" = ifelse(padj <= 0.05, TRUE, FALSE))

#30% VS Rest
deg.res <-
  des |>
  DESeq2::results(
    contrast = c(
      #DESeq2::resultsNames(des)
      "CO20.04" = - 1/2,
      "CO24.8" = - 1/2,
      "CO230" = + 1
    ),
    tidy = TRUE
  ) |>
  as_tibble() |>
  rename(Geneid = row) |>
  mutate(test = 'Other->30% CO2') |>
  mutate("is.de" = ifelse(padj <= 0.05, TRUE, FALSE)) |>
  bind_rows(deg.res) |>
  mutate(test = case_when(test == "4.8->30% CO2" ~ "4+8->30% CO2", test == "0.04->4.8% CO2" ~ "0.04->4+8% CO2", .default = test)) |>
  select(Geneid, baseMean, test, log2FoldChange, pvalue, padj, is.de) |>
  left_join(annot)

write_tsv(deg.res, paste0(output.dir, 'analysis/D_DEGs.tsv'))

#venn
deg.res |> mutate(log2FoldChange = round(log2FoldChange, 2)) |> filter(is.de & abs(log2FoldChange) >= 1) |> group_by(test) |> group_split() |> as.list() -> foo
set_names(
  foo |> map(pull, Geneid),
  foo |> map(pull, test) |> map(dplyr::first) |> unlist()) -> de.list
de.list |>
  venn::venn(
    ilabels = 'counts',
    ilcs = 1.6, sncs = 1.6,
    zcolor = cbPalette,
    box = FALSE,
    ggplot = TRUE
  )

deg30avg <- deg.res %>% filter(test == "Other->30% CO2")


################################################################################
# Normalized expression levels

des.vst <-
  des|>
  DESeq2::vst()

dat.vst <-
  des.vst |>
  SummarizedExperiment::assay()

dat.norm <-
  des |>
  DESeq2::counts(normalized = TRUE)

dat.vst |>
  as_tibble(rownames = 'Geneid') |>
  write_tsv(paste0(output.dir, 'analysis/D_vst-expression.tsv'))

dat.norm |>
  as_tibble(rownames = 'Geneid') |>
  write_tsv(paste0(output.dir, 'analysis/D_normalized-counts.tsv'))


################################################################################
# z-scaled gene expression

z.mat <-
  dat.vst |>
  apply(1, scale) |>
  t() |> 
  magrittr::set_colnames(colnames(dat.vst))

z.expr <-
  z.mat |>
  as_tibble(rownames = 'Geneid') |>
  pivot_longer(- Geneid, names_to = 'lib', values_to = 'z.expression') |>
  left_join(meta, 'lib') |>
  group_by(Geneid, CO2) |>
  summarize(avg = mean(z.expression)) |>
  ungroup() |>
  pivot_wider(names_from = 'CO2', values_from = 'avg')

z.expr |>
  write_tsv(paste0(output.dir, 'analysis/D_z-expr-mean-cond.tsv'))

z.expr |>
  rename('0.04% CO2' = '0.04',
         '4% CO2' = '4',
         '8% CO2' = '8',
         '30% CO2' = '30') |>
  inner_join(annot, by = "Geneid") |>
  transmute(Geneid, "0.04" = `0.04% CO2`, "4" =  `4% CO2`, "8" = `8% CO2`, "30" = `30% CO2`, "locus" = old_locus_tag) |>
  drop_na(locus) |>
  pivot_longer(cols = 2:5, names_to = "condition", values_to = "avg.z") |>
  readr::write_tsv(paste0(output.dir,"analysis/K_string-z-expression.tsv"))


################################################################################
# PCA

p.pca <- helper.pca(des.vst, sh = 'lane', color_values = cbPalette[c(1, 6, 2, 7)])
p.skree <- helper.skree(des.vst, nudge = 10) + xlim(c(1, 16))


p.pca | p.skree
ggsave(paste0(output.dir, 'analysis/D_PCA.jpeg'), width = 15, height = 8)


################################################################################
# Volcano plot

deg.res |>
  mutate(
    is.de = replace_na(is.de, FALSE),
    is.de = ifelse(is.de, 'yes' , 'no')
  ) |>
  ggplot(aes(log2FoldChange, -log10(padj), color = is.de)) +
  geom_point(alpha = .7, size = 1) +
  scale_color_manual(name = 'Diff. expressed', values = cbPalette[c(1, 7)]) +
  ylab('-log10 FDR') +
  facet_wrap(~ test) +
  theme_pubr(18)

ggsave(paste0(output.dir, 'analysis/D_volcano.jpeg'), width = 12, height = 8)


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
  breaks = seq(-2, 2, length.out = 60),
  color = colorRampPalette(rev(
    RColorBrewer::brewer.pal(n = 7, name = "RdBu")))(59),
  filename = paste0(output.dir, 'analysis/D_heatmap.jpeg')
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

jpeg(paste0(output.dir, 'analysis/D_logFC-cor.jpeg'), res = 400,
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
  write_tsv(paste0(output.dir, 'analysis/D_overview.tsv'))

#The largest numbers of DEGs (1,026 DEGs with abslogFC>=1) were observed for the comparisons of air to any other CO2 condition
dat.de |> mutate(log2FoldChange = round(log2FoldChange, 2)) |> filter((test == "0.04->30% CO2" | test == "0.04->4+8% CO2") & padj<=0.05 & abs(log2FoldChange)>=1) |> distinct(Geneid) |> dplyr::count()
#259 genes (abslogFC>=1) were differentially expressed in 30% compared to optimal growth condition
dat.de |> mutate(log2FoldChange = round(log2FoldChange, 2)) |> filter((test == "4+8->30% CO2") & padj<=0.05 & abs(log2FoldChange)>=1) |> distinct(Geneid) |> dplyr::count()
#732 genes (abslogFC>=1) were differentially expressed in 30% compared to at least one other condition (0.04% or 4/8% CO2) 
dat.de |> mutate(log2FoldChange = round(log2FoldChange, 2)) |> filter((test == "4+8->30% CO2" | test == "0.04->30% CO2") & padj<=0.05 & abs(log2FoldChange)>=1) |> distinct(Geneid) |> dplyr::count()
#220 genes (abslogFC>=1) were differentially expressed in 30% compared to the average of other conditions
dat.de |> mutate(log2FoldChange = round(log2FoldChange, 2)) |> filter((test == "Other->30% CO2") & padj<=0.05 & abs(log2FoldChange)>=1) |> distinct(Geneid) |> dplyr::count()


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
  write_tsv(paste0(output.dir, 'analysis/D_overview-by-type.tsv'))


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
      mutate(log2FoldChange = round(log2FoldChange, 2)) |>
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
  write_tsv(paste0(output.dir, 'analysis/D_overview-by-logFC.tsv'))


################################################################################
# Correlation between different contrasts

deg.res |>
  filter(test == "0.04->30% CO2" | test == "0.04->4+8% CO2" | test == "4+8->30% CO2") |>
  select(Geneid, test, log2FoldChange.pairwise = log2FoldChange) |>
  left_join(deg30avg, 'Geneid') |>
  ggscatter(
    'log2FoldChange', 'log2FoldChange.pairwise',
    add = 'reg.line', add.params = list(color = 'red'),
    cor.coef = TRUE, cor.coeff.args = list(color = 'red', size = 5),
    facet.by = 'test.x',
    alpha = .5
  ) +
  xlab('Other -> 30% CO2') +
  ylab('Pairwise comparison') +
  theme_pubr(18) +
  ggtitle('log2 Fold-Changes')


ggsave(paste0(output.dir, 'analysis/D_logFC-contrast-comparisons.jpeg'),
       width = 8, height = 8, dpi = 400)


################################################################################
# Heatmap of expression of DEGs Other->30%

deg30avg %>%
  filter(padj <= 0.05) |>
  mutate(log2FoldChange = round(log2FoldChange, 2)) |> 
  filter(abs(log2FoldChange) >= 1) |>
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
  border_color = NA,
  scale = 'none',
  cluster_cols = lib.clust,
  show_rownames = FALSE,
  show_colnames = FALSE,
  annotation_col = cl.df,
  annotation_colors = cl,
  color = colorRampPalette(rev(
    RColorBrewer::brewer.pal(n = 7, name = "RdBu")))(59),
  filename = paste0(output.dir, 'analysis/D_heatmap-30VsOther.jpeg')
)
dev.off()

# load as grid
p.heat <-
  ggdraw() +
  draw_image(
    paste0(output.dir, 'analysis/D_heatmap-30VsOther.jpeg')
  ) +
  theme_pubr(18) +
  theme(axis.line = element_blank())


################################################################################
# Compute gene set enrichment of KEGG pathways of DEGs Other->30% CO2

gsea <-
  deg30avg |>
  drop_na(old_locus_tag, log2FoldChange) |>
  arrange(desc(log2FoldChange)) |>
  with(set_names(log2FoldChange, old_locus_tag)) |>
  gseKEGG(
    organism = 'syp',
    maxGSSize = 70, # cutoff from comments above
    minGSSize = 5,
    pvalueCutoff = 0.05,
    eps = 0,
    seed = TRUE
  )
gsea@result$Description <- str_remove(gsea$Description, ' - Picosynechococcus sp. PCC 7002')

gsea |>
  as_tibble() |>
  write_tsv(paste0(output.dir, 'analysis/D_gsea-30VsOther.tsv'))

i <- gsea
p <- gseaplot2(i, 1:nrow(i), 
               base_size = 14,
               color = RColorBrewer::brewer.pal(nrow(i), 'Paired'))

# Change axis text for clarity
p[[3]] <-  p[[3]] + ylab('log2 Fold Change')
p
ggsave(paste0(output.dir, 'analysis/D_gsea-30VsOther.jpeg'),
       width = 8, height = 8, dpi = 400)


################################################################################
# input for network analysis in Cytoscape are DEGs Other->30%

annot |>
  filter(Geneid %in% mask) |>
  drop_na(old_locus_tag) |>
  pull(old_locus_tag) |>
  write_lines(paste0(output.dir, 'analysis/K_string-30VsOther.txt'))

deg30avg %>%
  filter(Geneid %in% mask) |>
  arrange(desc(log2FoldChange)) |>
  transmute(name, locus_tag, old_locus_tag, log2FoldChange, padj, product) |>
  write_tsv(paste0(output.dir, 'analysis/D_deg-30VsOther.tsv'))


################################################################################
# Overrepresentation analysis of KEGG pathways of DEGs Other->30% CO2

deg30avg <- deg30avg %>% mutate(Direction = if_else(`log2FoldChange` >= 0, true = "up", false = "down"))
deg30avg %>% mutate(log2FoldChange = round(log2FoldChange, 2)) |> filter(abs(log2FoldChange) >= 1 & padj <= 0.05) %>% select(Direction) %>% table()

oa.30VsOther.up <- enrichKEGG(gene = deg30avg |> filter(Direction == "up") |> select(old_locus_tag) |> unlist(),
           organism = 'syp',
           pAdjustMethod = "BH", qvalueCutoff  = 0.05)
oa.30VsOther.down <- enrichKEGG(gene = deg30avg |> filter(Direction == "down") |> select(old_locus_tag) |> unlist(),
                               organism = 'syp',
                               pAdjustMethod = "BH", qvalueCutoff  = 0.05)
oa.30VsOther.up@result$direction <- "Upregulated"
oa.30VsOther.down@result$direction <- "Downregulated"
combined <- rbind(oa.30VsOther.up@result, oa.30VsOther.down@result)
ggplot(combined[combined$p.adjust<0.05,], aes(x = direction, y = Description)) +
  geom_point(aes(size = Count, color = direction)) +
  scale_color_manual(values = c("Upregulated" = "red", "Downregulated" = "blue")) +
  theme_minimal() +
  labs(title = "Significant KEGG Pathways",
       x = "Regulation Direction",
       y = "Pathway",
       color = "Direction")
ggsave(paste0(output.dir, 'analysis/D_oa-30VsOther.jpeg'),
      width = 8, height = 8, dpi = 400)
dotplot(oa.30VsOther.down, showCategory = 10, label_format = 50) + 
  theme(axis.text.y = element_text(size = 8))
ggsave(paste0(output.dir, 'analysis/D_oa-30VsOther-dotplot.jpeg'),
       width = 8, height = 8, dpi = 400)


################################################################################
# Compute gene set enrichment of KEGG pathways for 30% vs optimal growth (4+8% CO2)

t.co2 <-
  deg.res |>
  filter(test == "4+8->30% CO2" & !is.na(locus_tag)) |>
  mutate(old_locus_tag = ifelse(is.na(old_locus_tag), locus_tag, old_locus_tag)) |>
  select('old_locus_tag', 'log2FoldChange', 'padj', 'is.de') |>
  rename('Locus tag' = 'old_locus_tag', 'log2FoldChange-HighCO2' = 'log2FoldChange')
t.co2 <- t.co2 %>% mutate(Direction = if_else(`log2FoldChange-HighCO2` >= 0, true = "up", false = "down"))
t.co2 %>% dplyr::filter(abs(`log2FoldChange-HighCO2`) > 1 & padj <= 0.05) |> select(Direction) |> table()
gsea <-
  t.co2 |>
  drop_na(`Locus tag`, `log2FoldChange-HighCO2`) |>
  arrange(desc(`log2FoldChange-HighCO2`)) |>
  with(set_names(`log2FoldChange-HighCO2`, `Locus tag`)) |>
  gseKEGG(
    organism = 'syp',
    maxGSSize = 70, # cutoff from comments above
    minGSSize = 5,
    pvalueCutoff = 0.05,
    eps = 0,
    seed = TRUE
  )
gsea@result$Description <- str_remove(gsea$Description, ' - Picosynechococcus sp. PCC 7002')

gsea |>
  as_tibble() |>
  write_tsv(paste0(output.dir, 'analysis/D_gsea-30VsOpt.tsv'))

i <- gsea
p <- gseaplot2(i, 1:nrow(i), 
               base_size = 14,
               color = RColorBrewer::brewer.pal(nrow(i), 'Paired'))
# Change axis text for clarity
p[[3]] <-  p[[3]] + ylab('log2 Fold Change')
p
ggsave(paste0(output.dir, 'analysis/D_gsea-30VsOpt.jpeg'),
       width = 8, height = 8, dpi = 400)

gsea@result$core_enrichment[gsea@result$Description=="Nitrogen metabolism"]
gsea@result$core_enrichment[gsea@result$Description=="Arginine biosynthesis"]


################################################################################
# Compute gene set enrichment of KEGG pathways for 30% vs 0.04% CO2

t.co2 <-
  deg.res |>
  filter(test == "0.04->30% CO2" & !is.na(locus_tag)) |>
  mutate(old_locus_tag = ifelse(is.na(old_locus_tag), locus_tag, old_locus_tag)) |>
  select('old_locus_tag', 'log2FoldChange', 'padj', 'is.de') |>
  rename('Locus tag' = 'old_locus_tag', 'log2FoldChange-HighCO2' = 'log2FoldChange')
t.co2 <- t.co2 %>% mutate(Direction = if_else(`log2FoldChange-HighCO2` >= 0, true = "up", false = "down"))
t.co2 %>% dplyr::filter(abs(`log2FoldChange-HighCO2`) > 1 & padj <= 0.05) |> select(Direction) |> table()
gsea <-
  t.co2 |>
  drop_na(`Locus tag`, `log2FoldChange-HighCO2`) |>
  arrange(desc(`log2FoldChange-HighCO2`)) |>
  with(set_names(`log2FoldChange-HighCO2`, `Locus tag`)) |>
  gseKEGG(
    organism = 'syp',
    maxGSSize = 70, # cutoff from comments above
    minGSSize = 5,
    pvalueCutoff = 0.05,
    eps = 0,
    seed = TRUE
  )
gsea@result$Description <- str_remove(gsea$Description, ' - Picosynechococcus sp. PCC 7002')

gsea |>
  as_tibble() |>
  write_tsv(paste0(output.dir, 'analysis/D_gsea-30VsAir.tsv'))

i <- gsea
p <- gseaplot2(i, 1:nrow(i), 
               base_size = 14,
               color = RColorBrewer::brewer.pal(nrow(i), 'Paired'))
# Change axis text for clarity
p[[3]] <-  p[[3]] + ylab('log2 Fold Change')
p
ggsave(paste0(output.dir, 'analysis/D_gsea-30VsAir.jpeg'),
       width = 8, height = 8, dpi = 400)
gsea@result$core_enrichment[gsea@result$Description=="Sulfur metabolism"]
gsea@result$core_enrichment[gsea@result$Description=="Folate biosynthesis"]


################################################################################
# Compute gene set enrichment of KEGG pathways for 8+4% vs 0.04% CO2

t.co2 <-
  deg.res |>
  filter(test == "0.04->4+8% CO2" & !is.na(locus_tag)) |>
  mutate(old_locus_tag = ifelse(is.na(old_locus_tag), locus_tag, old_locus_tag)) |>
  select('old_locus_tag', 'log2FoldChange', 'padj', 'is.de') |>
  rename('Locus tag' = 'old_locus_tag', 'log2FoldChange-HighCO2' = 'log2FoldChange')
t.co2 <- t.co2 %>% mutate(Direction = if_else(`log2FoldChange-HighCO2` >= 0, true = "up", false = "down"))
t.co2 %>% dplyr::filter(abs(`log2FoldChange-HighCO2`) > 1 & padj <= 0.05) |> select(Direction) |> table()
gsea <-
  t.co2 |>
  drop_na(`Locus tag`, `log2FoldChange-HighCO2`) |>
  arrange(desc(`log2FoldChange-HighCO2`)) |>
  with(set_names(`log2FoldChange-HighCO2`, `Locus tag`)) |>
  gseKEGG(
    organism = 'syp',
    maxGSSize = 70, # cutoff from comments above
    minGSSize = 5,
    pvalueCutoff = 0.05,
    eps = 0,
    seed = TRUE
  )
gsea@result$Description <- str_remove(gsea$Description, ' - Picosynechococcus sp. PCC 7002')

gsea |>
  as_tibble() |>
  write_tsv(paste0(output.dir, 'analysis/D_gsea-OptVsAir.tsv'))

i <- gsea
p <- gseaplot2(i, 1:nrow(i), 
               base_size = 14,
               color = RColorBrewer::brewer.pal(nrow(i), 'Paired'))
# Change axis text for clarity
p[[3]] <-  p[[3]] + ylab('log2 Fold Change')
p
ggsave(paste0(output.dir, 'analysis/D_gsea-OptVsAir.jpeg'),
       width = 8, height = 8, dpi = 400)
ggsave(paste0(output.dir, 'analysis/D_gsea-OptVsAir.svg'),
       width = 7, height = 10, dpi = 400)


################################################################################
# Overrepresentation analysis of KEGG pathways of DEGs from pairwise comparisons and their logical combinations

deg.30VsAir <- deg.res |>
  filter(test == "0.04->30% CO2") |>
  mutate(Direction = if_else(`log2FoldChange` >= 0, true = "up", false = "down")) |>
  mutate(log2FoldChange = round(log2FoldChange, 2)) |>
  filter(abs(log2FoldChange) >= 1 & padj <= 0.05)

deg.30VsOpt <- deg.res |>
  filter(test == "4+8->30% CO2") |>
  mutate(Direction = if_else(`log2FoldChange` >= 0, true = "up", false = "down")) |>
  mutate(log2FoldChange = round(log2FoldChange, 2)) |>
  filter(abs(log2FoldChange) >= 1 & padj <= 0.05)

deg.OptVsAir <- deg.res |>
  filter(test == "0.04->4+8% CO2") |>
  mutate(Direction = if_else(`log2FoldChange` >= 0, true = "up", false = "down")) |>
  mutate(log2FoldChange = round(log2FoldChange, 2)) |>
  filter(abs(log2FoldChange) >= 1 & padj <= 0.05)

# DEGs 30% vs 0.04% CO2
oa.30VsAir.up <- enrichKEGG(gene = deg.30VsAir |> filter(Direction == "up") |> select(old_locus_tag) |> unlist(),
                            organism = 'syp',
                            pAdjustMethod = "BH", qvalueCutoff  = 0.05)
oa.30VsAir.down <- enrichKEGG(gene = deg.30VsAir |> filter(Direction == "down") |> select(old_locus_tag) |> unlist(),
                            organism = 'syp',
                            pAdjustMethod = "BH", qvalueCutoff  = 0.05)
oa.30VsAir.up@result$direction <- "Upregulated"
oa.30VsAir.down@result$direction <- "Downregulated"
combined <- rbind(oa.30VsAir.up@result, oa.30VsAir.down@result)
ggplot(combined[combined$p.adjust<0.05,], aes(x = direction, y = Description)) +
  geom_point(aes(size = Count, color = direction)) +
  scale_color_manual(values = c("Upregulated" = "red", "Downregulated" = "blue")) +
  theme_minimal() +
  labs(title = "Significant KEGG Pathways",
       x = "Regulation Direction",
       y = "Pathway",
       color = "Direction")
ggsave(paste0(output.dir, 'analysis/D_oa-30VsAir.jpeg'),
       width = 8, height = 8, dpi = 400)

# DEGs 30% vs optimal growth
oa.30VsOpt.up <- enrichKEGG(gene = deg.30VsOpt |> filter(Direction == "up") |> select(old_locus_tag) |> unlist(),
                            organism = 'syp',
                            pAdjustMethod = "BH", qvalueCutoff  = 0.05)
oa.30VsOpt.down <- enrichKEGG(gene = deg.30VsOpt |> filter(Direction == "down") |> select(old_locus_tag) |> unlist(),
                              organism = 'syp',
                              pAdjustMethod = "BH", qvalueCutoff  = 0.05)
oa.30VsOpt.up@result$direction <- "Upregulated"
oa.30VsOpt.down@result$direction <- "Downregulated"
combined <- rbind(oa.30VsOpt.up@result, oa.30VsOpt.down@result)
ggplot(combined[combined$p.adjust<0.05,], aes(x = direction, y = Description)) +
  geom_point(aes(size = Count, color = direction)) +
  scale_color_manual(values = c("Upregulated" = "red", "Downregulated" = "blue")) +
  theme_minimal() +
  labs(title = "Significant KEGG Pathways",
       x = "Regulation Direction",
       y = "Pathway",
       color = "Direction")
ggsave(paste0(output.dir, 'analysis/D_oa-30VsOpt.jpeg'),
       width = 8, height = 8, dpi = 400)

# DEGs 8+4% vs 0.04% CO2
oa.OptVsAir.up <- enrichKEGG(gene = deg.OptVsAir |> filter(Direction == "up") |> select(old_locus_tag) |> unlist(),
                            organism = 'syp',
                            pAdjustMethod = "BH", qvalueCutoff  = 0.05)
oa.OptVsAir.down <- enrichKEGG(gene = deg.OptVsAir |> filter(Direction == "down") |> select(old_locus_tag) |> unlist(),
                              organism = 'syp',
                              pAdjustMethod = "BH", qvalueCutoff  = 0.05)
oa.OptVsAir.up@result$direction <- "Upregulated"
oa.OptVsAir.down@result$direction <- "Downregulated"
combined <- rbind(oa.OptVsAir.up@result, oa.OptVsAir.down@result)
ggplot(combined[combined$p.adjust<0.05,], aes(x = direction, y = Description)) +
  geom_point(aes(size = Count, color = direction)) +
  scale_color_manual(values = c("Upregulated" = "red", "Downregulated" = "blue")) +
  theme_minimal() +
  labs(title = "Significant KEGG Pathways",
       x = "Regulation Direction",
       y = "Pathway",
       color = "Direction")
ggsave(paste0(output.dir, 'analysis/D_oa-OptVsAir.jpeg'),
       width = 8, height = 8, dpi = 400)

# DEGs (30% vs 0.04% CO2) OR (30% vs optimal growth)
oa.30VsAirOR30VsOpt <- enrichKEGG(gene = deg.30VsAir |> filter(!is.na(old_locus_tag)) |> select(old_locus_tag) |> full_join(deg.30VsOpt |> filter(!is.na(old_locus_tag)) |> select(old_locus_tag), by = "old_locus_tag") |> unlist(),
                            organism = 'syp',
                            pAdjustMethod = "BH", qvalueCutoff  = 0.05)
dotplot(oa.30VsAirOR30VsOpt)
ggsave(paste0(output.dir, 'analysis/D_oa-30VsAirOR30VsOpt.jpeg'),
       width = 8, height = 8, dpi = 400)

# DEGs (30% vs 0.04% CO2) AND (30% vs optimal growth)
oa.30VsAirAND30VsOpt <- enrichKEGG(gene = deg.30VsAir |> filter(!is.na(old_locus_tag)) |> select(old_locus_tag) |>
                                     inner_join(deg.30VsOpt |> filter(!is.na(old_locus_tag)) |> select(old_locus_tag), by = "old_locus_tag") |> unlist(),
                                  organism = 'syp',
                                  pAdjustMethod = "BH", qvalueCutoff  = 0.05)
dotplot(oa.30VsAirAND30VsOpt)
ggsave(paste0(output.dir, 'analysis/D_oa-30VsAirAND30VsOpt.jpeg'),
       width = 8, height = 8, dpi = 400)

# DEGs (30% vs 0.04% CO2) AND (30% vs optimal growth) AND NOT (8+4% vs 0.04% CO2)
oa.30VsAirAND30VsOptANDNOTOptVsAir <- enrichKEGG(gene = deg.30VsAir |> filter(!is.na(old_locus_tag)) |> select(old_locus_tag) |>
                                                   inner_join(deg.30VsOpt |> filter(!is.na(old_locus_tag)) |> select(old_locus_tag), by = "old_locus_tag") |>
                                                   anti_join(deg.OptVsAir |> filter(!is.na(old_locus_tag)) |> select(old_locus_tag), by = "old_locus_tag") |>
                                                   unlist(),
                                   organism = 'syp',
                                   pAdjustMethod = "BH", qvalueCutoff  = 0.05)
dotplot(oa.30VsAirAND30VsOptANDNOTOptVsAir)
ggsave(paste0(output.dir, 'analysis/D_oa-30VsAirAND30VsOptANDNOTOptVsAir.jpeg'),
       width = 8, height = 8, dpi = 400)

# The "other vs 30%" comparison misses genes that get strongly down-regulated in one pairwise comparison of 0.04% or optimal growth to 30% and
# strongly up-regulated in the other pairwise comparison to 30%:
# 38 genes
deg.30VsOther <- deg.res |>
  filter(test == "Other->30% CO2") |>
  mutate(Direction = if_else(`log2FoldChange` >= 0, true = "up", false = "down")) |>
  mutate(log2FoldChange = round(log2FoldChange, 2)) |>
  filter(abs(log2FoldChange) >= 1 & padj <= 0.05)

deg.30VsAir |> inner_join(deg.30VsOpt, by = "Geneid") |>
  anti_join(deg.30VsOther, by = "Geneid") |>
  filter(Direction.x != Direction.y) |>
  write_tsv(paste0(output.dir, 'analysis/D_deg-30VsOptAND30VsAirANDNOT30VsOther.tsv'))


################################################################################
################################################################################
sessionInfo()
