# Investigate the impact of differential expression on amino acid compositions

set.seed(123)


library(tidyverse)
library(ggpubr)
library(patchwork)

library(plyranges)
library(Biostrings)
library(BSgenome)

library(conflicted)
conflicts_prefer(plyranges::filter)
conflicts_prefer(dplyr::rename)
conflicts_prefer(base::intersect)
conflicts_prefer(dplyr::n)

conflicts_prefer(purrr::map)
conflicts_prefer(plyranges::select)

# nicer amino acid names
data(aaMap, package = 'Biobase')

# https://riptutorial.com/r/example/28354/colorblind-friendly-palettes
cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442",
               "#0072B2", "#D55E00", "#CC79A7")

################################################################################
################################################################################
# Load data

annot <-
  'data/C_annotation.tsv' |>
  read_tsv()

coord <-
  'raw-data/co2-adaptation-counts.tsv' |>
  read_tsv() |>
  select(Geneid, Chr, Start, End, Strand, Length)

deg <-
  'analysis/D_stagewise-adjusted-DEGs.tsv' |>
  read_tsv()

norm <-
  'analysis/D_normalized-counts.tsv' |>
  read_tsv()

vst <-
  'analysis/D_vst-expression.tsv' |>
  read_tsv()

meta <-
  'data/C_meta.tsv' |>
  read_tsv()

genome <-
  'raw-data/PCC7002-genome.fna.gz' |>
  readDNAStringSet()
# only accession number in name
names(genome) <- genome |> names() |> str_remove(' .*$')


################################################################################
# Sequences of coding genes and amino acid composition

cds.range <-
  annot |>
  filter(type == 'protein_coding') |>
  select(Geneid, product, name) |>
  left_join(coord, 'Geneid') |>
  select(- Length) |>
  rename(
    seqnames = Chr,
    start = Start, end = End,
    strand = Strand
  ) |>
  as_granges()

# extract sequence and convert to amino acids
nts <- getSeq(genome, cds.range)
bcode <- getGeneticCode('Bacterial, Archaeal and Plant Plastid',
                        full.search = TRUE)
aas <- translate(nts, genetic.code = bcode)

# find 'complete' peptides
peptide.mask <- 
  # starts with methionine
  (subseq(aas, 1, 1) == 'M') &
  # only peptides, no odd undefined AAs
  (subseq(aas, 2, -2) %>% as.character() %>% str_detect('^[A-Z]+$')) &
  # ends with termination
  (subseq(aas, -1, -1) == '*')
# table(peptide.mask) / length(cds.range) * 100
# FALSE        TRUE 
# 0.03137747 99.96862253 

# Only complete peptides, excl.start and stop codon
aas2 <- aas[peptide.mask]
names(aas2) <- cds.range$Geneid[peptide.mask]
aas2 <- Biostrings::subseq(aas2, 2, -2)

writeXStringSet(aas2, 'analysis/G_peptides.faa')

# The frequencies
freqs <- alphabetFrequency(aas2, as.prob = TRUE)#* 100
rownames(freqs) <- names(aas2)

# exclude all 0s, and keep only DEGs
mask1 <- deg |>
  filter(is.de) |>
  pull(Geneid) |>
  intersect(rownames(freqs))
mask2 <- apply(freqs, 2, max) > 0
freqs2 <- freqs[mask1, mask2]
# nicer AA names
colnames(freqs2) <- with(aaMap, set_names(name, let.1))[colnames(freqs2)]

# corresponding vst expression value and gene lengths
vst2 <- vst.mat[rownames(freqs2), ]
len2 <- with(de.len, set_names(Length, Geneid))[rownames(freqs2)] |> log10()

# export freqs
freqs2 |>
  as_tibble(rownames = 'Geneid') |>
  write_tsv('analysis/G_frequencies.tsv')

################################################################################
# Overall AA freq

freqs2 |>
  as_tibble(rownames = 'Geneid') |>
  pivot_longer(- Geneid) |>
  mutate(AA = fct_reorder(name, value)) |>
  ggplot(aes(AA, value, fill = AA)) +
  geom_violin() +
  # geom_boxplot(fill = 'white', width = .2) +
  xlab(NULL) +
  ylab('Gene length normalized\nAmino acid frequency') +
  theme_pubr(18) +
  theme(
    legend.position = 'hide',
    axis.text.x = element_text(angle = 60, hjust = 1)
  )

ggsave('analysis/G_freqs-overall.jpeg', width = 10, height = 7, dpi = 400)

################################################################################
# Explore distribution of AA frequencies and potential connection to gene lengths

# Z-scale frequency per amino acid
f2z <-
  freqs2 |>
  sqrt() |>
  apply(2, scale) |>
  magrittr::set_rownames(rownames(freqs2))

# Check for "normalnesss"
list(
  'AA freq.' = freqs2,
  'sqrt AA freq.' = sqrt(freqs2),
  'z-scaled sqrt AA freq.' = f2z
) %>%
  map2(names(.), ~ tibble(x = c(.x), n = .y)) |>
  bind_rows() |>
  mutate_at('n', fct_inorder) -> foo

p3.hist <-
  foo |>
    gghistogram('x', facet.by = 'n', scales = 'free') +
    ggtitle('Histograms') +
    theme_pubr(18)
p3.qq <-
  foo |>
    group_by(n) |>
    reframe(
      qs = seq(0, 1, .001),
      qo = quantile(x, probs = qs),
      qn = qnorm(qs)
    ) |>
    ggscatter('qn', 'qo', facet.by = 'n', scales = 'free') +
    geom_abline(slope = 1, color = 'blue') +
    xlab('Normal distribution') +
    ylab('Observed values') +
    ggtitle('Quantile-quantile plot') +
    theme_pubr(18)

# Correlation length ~ frequency
p3.cor <-
  f2z |>
  as_tibble(rownames = 'Geneid') |>
  left_join(coord |> select(Geneid, Length), 'Geneid') |>
  mutate_at('Length', log10) |>
  pivot_longer(- c(Geneid, Length)) |>
  ggplot(aes(Length, value)) +
  geom_hex() +
  scale_fill_viridis_c() +
  guides(fill = guide_colorbar(barwidth = unit(5, 'cm'))) +
  geom_smooth(method = 'lm', color = 'red', se = FALSE) +
  stat_cor(color = 'red', size = 5, label.x.npc = "middle") +
  scale_x_log10() +
  xlab('Gene length, log') +
  ylab('z-scaled sqrt AA freq.') +
  theme_pubr(18)


p3 <-
  ( p3.hist / p3.qq / p3.cor ) +
  plot_layout(heights = c(1, 1, 2)) +
  plot_annotation(tag_levels = 'A')

ggsave('analysis/G_AA-dist-cor.jpeg', p3,
       width = 12, height = 12, dpi = 400)

################################################################################
# Correlate AA and expression

# Per gene z-Scaled expression
vz <-
  vst2 |>
  apply(1, scale) |>
  t() |>
  magrittr::set_colnames(colnames(vst2))

# # keep only top X most variant expressed genes among the DEGs
# v.dat <-
#   tibble(
#     gene = rownames(vst2),
#     v =  MatrixGenerics::rowVars(vst2)
#   )
# 
# qv <- quantile(v.dat$v, .9)
# nice <-
#   sprintf(
#     'Top 10%% (n = %g genes) with\nvariance ≥ %.2f',
#     v.dat |>
#       filter(v >= qv) |>
#       nrow(),
#     qv
#   )
# 
# v.dat |>
#   ggplot(aes( v)) +
#   stat_ecdf() +
#   geom_hline(yintercept = .9, linetype = 'dotted') +
#   geom_vline(xintercept = qv, color = 'red') +
#   scale_y_continuous(breaks = seq(0, 1, .1)) +
#   annotate(
#     'text', label = nice,
#     x = 1, y = .84, color = 'red',
#     size = 8,
#     hjust = 0
#   ) +
#   xlab('Expression variance (vst)') +
#   ylab('Empirical cumulative density') +
#   theme_pubr(18)
# 
# ggsave('analysis/I_top-var.jpeg',
#        width = 8, height = 8, dpi = 400)
# 
# 
# top.var <-
#   v.dat |>
#   filter(v >= qv) |>
#   pull(gene)

# nicer library names
look <-
  meta |>
  mutate(nice = sprintf('%s%% (%s)', CO2, sample)) |>
  with(set_names(nice, lib))

deg |>
  filter(is.de, abs(log2FoldChange) >= 1) |>
  pull(Geneid) |>
  unique() |>
  intersect(rownames(f2z)) -> mask

crossing(
  AA = colnames(f2z),
  lib = colnames(vst2)
) |>
  group_by(AA, lib) |>
  # do(r2 = cor.test(f2z[top.var, .$AA], vz[top.var, .$lib])$estimate) |>
  # do(r2 = cor.test(f2z[, .$AA], vz[, .$lib])$estimate) |>
  do(r2 = cor.test(f2z[mask, .$AA], vst2[mask, .$lib])$estimate) |>
  # do(r2 = cor.test(sqrt(freqs2)[mask, .$AA], vz[mask, .$lib])$estimate) |>
  ungroup() |>
  unnest(r2) |>
  mutate(lib = look[lib]) |>
  spread(lib, r2) -> foo

foo |>
  select(- AA) |>
  as.matrix() |>
  magrittr::set_rownames(foo$AA) |>
  pheatmap::pheatmap(
    scale = 'none',
    show_colnames = TRUE,
    show_rownames = TRUE,
    display_numbers = TRUE,
    fontsize = 18,
    number_color = 'black',
    color = colorRampPalette(rev(
      RColorBrewer::brewer.pal(n = 5, name = "RdBu")))(59),
    # filename = 'analysis/I_AA-expr-cor.jpeg', width = 14, height = 14
  )
dev.off()
  


################################################################################
# Signal in frequencies

pca <- prcomp(t(f2z))
percentVar <- ( pca$sdev^2/sum(pca$sdev^2) ) %>%
  set_names(paste0('PC', 1:length(.)))

x <- 'PC1'
y <- 'PC2'
pca$x |>
  as_tibble(rownames = 'amino') |>
  left_join(aaMap, c('amino' = 'name')) |>
  ggplot(aes(!! sym(x), !! sym(y), color = scProp)) +
  # ggplot(aes(!! sym(x), !! sym(y))) +
  geom_point(size = 5) +
  ggsci::scale_color_jama(name = 'side chain property at pH 7') +
  xlab(sprintf('%s: %.1f%% variance', x, percentVar[[x]] * 100)) +
  ylab(sprintf('%s: %.1f%% variance', y, percentVar[[y]] * 100)) +
  ggrepel::geom_label_repel(aes(label = amino), show.legend = FALSE, size = 5) +
  coord_cartesian() +
  theme_pubr(18) +
  ggtitle('PCA of amino acid frequencies alone')

ggsave('analysis/I_AA-freq-pca.jpeg', width = 10, height = 8, dpi = 400)

################################################################################
################################################################################
