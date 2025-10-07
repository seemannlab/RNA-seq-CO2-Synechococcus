library(DESeq2)
library(tidyverse)
library(ggplot2)
library(readxl)
library(pheatmap)
library(venn)
library(ggrepel)
library(grid)
library(ggpubr)

# https://riptutorial.com/r/example/28354/colorblind-friendly-palettes
cbPalette <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442",
               "#0072B2", "#D55E00", "#CC79A7", "#999999")

analysis.dir <- '/media/dna/projects/rth/co2capture/continuation-ses/subprojects/RNA-seq_1v15-pct-CO2/RNA-seq-CO2-Synechococcus/stress/'
out.dir <- "/media/dna/projects/rth/co2capture/continuation-ses/subprojects/RNA-seq_1v15-pct-CO2/Manuscript-Overleaf/files-manuscript/"
supp.outdir <- "/media/dna/projects/rth/co2capture/continuation-ses/subprojects/RNA-seq_1v15-pct-CO2/Manuscript-Overleaf/files-supplement/"

#Ludwig and Bryant 2011
f.ludwigbryant2011 <- file.path(analysis.dir, 'LudwigBryant2011_Data Sheet 2.XLSX')
t.ludwigbryant2011 <- read_xlsx(f.ludwigbryant2011, range = cell_rows(4:3190))

#filter high light
t.ludwigbryant2011 <-
  t.ludwigbryant2011 %>% select(`Locus tag`,`Ratio high light/ standard`, `p value...3`) |>
  mutate(`p value...3` = ifelse(`p value...3`=="<0.0001", 0.0001, as.numeric(`p value...3`))) |>
  mutate(is.de = ifelse(`p value...3` <= 0.05, 'yes' , 'no')) |>
  mutate("log2FoldChange-HighLight" = ifelse(`Ratio high light/ standard`==0, 0, log2(as.numeric(`Ratio high light/ standard`))))
#check if all log2FC cutoffs are significant differential
t.ludwigbryant2011 %>% filter(`log2FoldChange-HighLight` >= 1) %>% select(is.de) %>% table()
#is.de
#no yes 
#1 224


#Ludwig and Bryant 2012A
f.ludwigbryant2012a <- file.path(analysis.dir, 'LudwigBryant2012A_Data Sheet 2.XLSX')
t.ludwigbryant2012a <- read_xlsx(f.ludwigbryant2012a, range = cell_rows(4:3190))

#filter heat and high salt
t.ludwigbryant2012a <-
  t.ludwigbryant2012a |>
  select(`Locus tag`, `Ratio heat shock/ standard (S4)`, `p value...15`, `Ratio 1.5M NaCl/ standard (S4)`, `p value...13`) |>
  mutate(`p value...15` = ifelse(`p value...15`=="<0.0001", 0.0001, as.numeric(`p value...15`))) |>
  mutate(is.de...15 = ifelse(`p value...15` <= 0.05, 'yes' , 'no')) |>
  mutate("log2FoldChange-Heat" = ifelse(`Ratio heat shock/ standard (S4)`==0, 0, log2(as.numeric(`Ratio heat shock/ standard (S4)`)))) |>
  mutate(`p value...13` = ifelse(`p value...13`=="<0.0001", 0.0001, as.numeric(`p value...13`))) |>
  mutate(is.de...13 = ifelse(`p value...13` <= 0.05, 'yes' , 'no')) |>
  mutate("log2FoldChange-HighSalt" = ifelse(`Ratio 1.5M NaCl/ standard (S4)`==0, 0, log2(as.numeric(`Ratio 1.5M NaCl/ standard (S4)`))))
#check if all log2FC cutoffs are significant differential
t.ludwigbryant2012a %>% filter(`log2FoldChange-Heat` >= 1) %>% select(is.de...15) %>% table()
#is.de...15
#no yes 
#8 359
t.ludwigbryant2012a %>% filter(`log2FoldChange-HighSalt` >= 1) %>% select(is.de...13) %>% table()
#is.de...13
#no yes 
#22 650

#Ludwig and Bryant 2012B
f.ludwigbryant2012b <- file.path(analysis.dir, 'LudwigBryant2012B_Data Sheet 3.XLSX')
t.ludwigbryant2012b <- read_xlsx(f.ludwigbryant2012b, range = cell_rows(4:3190))

#filter N limitation
t.ludwigbryant2012b <- 
  t.ludwigbryant2012b |>
  select(Locus_tag, `Ratio N-limited/ standard`, `p value...5`) |>
  mutate(`p value...5` = ifelse(`p value...5`=="<0.0001", 0.0001, as.numeric(`p value...5`))) |>
  mutate(is.de = ifelse(`p value...5` <= 0.05, 'yes' , 'no')) |>
  mutate("log2FoldChange-LowN" = ifelse(`Ratio N-limited/ standard`==0, 0, log2(as.numeric(`Ratio N-limited/ standard`)))) |>
  dplyr::rename('Locus tag' = 'Locus_tag')
#check if all log2FC cutoffs are significant differential
t.ludwigbryant2012b %>% filter(`log2FoldChange-LowN` >= 1) %>% select(is.de) %>% table()
#is.de
#no yes 
#2 590

#CO2 levels
f.co2 <- file.path(analysis.dir, '../analysis/D_DEGs.tsv')
t.co2 <- read.csv(f.co2, , sep = '\t', header = T)

#filter CO2
t.co2 <-
  t.co2 |>
  filter(test == "4+8->30% CO2" & !is.na(locus_tag)) |>
  mutate(old_locus_tag = ifelse(is.na(old_locus_tag), locus_tag, old_locus_tag)) |>
  select('old_locus_tag', 'log2FoldChange', 'padj', 'is.de') |>
  dplyr::rename('Locus tag' = 'old_locus_tag', 'log2FoldChange-HighCO2' = 'log2FoldChange')

#check if all log2FC cutoffs are significant differential
t.co2 %>% filter(`log2FoldChange-HighCO2` >= 1) %>% select(is.de) %>% table()
#is.de
#FALSE  TRUE 
#1   142

#join tables
t.stress <- 
  t.co2 |>
  inner_join(t.ludwigbryant2011, "Locus tag") |>
  inner_join(t.ludwigbryant2012a, "Locus tag") |>
  inner_join(t.ludwigbryant2012b, "Locus tag")

#filter genes with log2FC>1 in at least one stress condition
t.stress <- 
  t.stress |>
  #filter(abs(`log2FoldChange-HighCO2`)>1 | abs(`log2FoldChange-HighLight`)>1 | abs(`log2FoldChange-Heat`)>1 | abs(`log2FoldChange-HighSalt`)>1 | abs(`log2FoldChange-LowN`)>1) |>
  filter((`log2FoldChange-HighCO2`>=1 & is.de.x) | `log2FoldChange-HighLight`>=1 | `log2FoldChange-Heat`>=1 | `log2FoldChange-HighSalt`>=1 | `log2FoldChange-LowN`>=1) |>
  select(`Locus tag`, `log2FoldChange-HighCO2`, `log2FoldChange-HighLight`, `log2FoldChange-Heat`, `log2FoldChange-HighSalt`, `log2FoldChange-LowN`) |>
  column_to_rownames(var = "Locus tag") |>
  rename_with(~ str_replace(., "log2FoldChange-", ""))

#PCA
replace_non_finite <- function(df) {
  # Check if the input is a dataframe
  if (!is.data.frame(df)) {
    stop("Input must be a dataframe.")
  }
  
  # Apply the check and replacement to each element
  for (i in 1:ncol(df)) {
    if (is.numeric(df[, i])) { # Only process numeric columns
      for (j in 1:nrow(df)) {
        if (!is.finite(df[j, i])) {
          df[j, i] <- 0
        }
      }
    }
  }
  return(df)
}
t.stress <- replace_non_finite(t.stress)

pca <- prcomp(t(t.stress), scale. = TRUE)
pca_data <- data.frame(PC1 = pca$x[, 1], PC2 = pca$x[, 2], PC3 = pca$x[, 3], PC4 = pca$x[, 4])

# Calculate the proportion of variance explained:
variance_explained <- pca$sdev^2 / sum(pca$sdev^2)
scree_data <- data.frame(PC = 1:length(variance_explained), Variance = variance_explained)
stress.scree <- ggplot(scree_data, aes(x = PC, y = Variance)) +
  geom_line() +
  geom_point() +
  labs(title = "Scree Plot", x = "Principal Component", y = "Proportion") +
  theme_pubr(16)

#plot PCA
pca_data$labels <- colnames(t.stress)
pca_data$labels <- c("High~CO[2]", "High~Light", "Heat", "High~Salt", "Limited~N")
stress.pca <- ggplot(pca_data, aes(x = PC1, y = PC2, label = labels)) +
  geom_point() +
  geom_text(size=7, position=position_nudge(x = 2.6, y = 2), parse = TRUE) +
  labs(title = "PCA Plot", x = paste0("PC1 (", round(100*variance_explained[1],1), "%)"), y = paste0("PC2 (", round(100*variance_explained[2],1), "%)")) +
  theme_pubr(20)

#hierarchical clustering
clust.stress <- t.stress |> t() |> dist(method="euclidian") |> hclust()

#heat map
t.stress.sortedco2 <- t.stress[order(t.stress$HighCO2, decreasing = TRUE), ]
colnames(t.stress.sortedco2) <- c("High CO2", "High Light", "Heat", "High Salt", "Limited N")
stress.heatmap <- pheatmap(
  t.stress.sortedco2[1:50,],
  #scale = 'column',
  scale = 'none',
  show_rownames = TRUE,
  cluster_cols = clust.stress,
  cluster_rows = FALSE,
  breaks = seq(-1, 3, length.out = 10),
  color = colorRampPalette(rev(
    RColorBrewer::brewer.pal(n = 4, name = "RdBu")))(9),
  angle_col = 45,
  fontsize_col = 18,
  silent = TRUE)

#save log2fc of stress proteins
readr::write_tsv(t.stress.sortedco2 |> rownames_to_column(var = 'Geneid'), file.path(analysis.dir, 'tbl.stress.tsv'))

#jpeg(paste0(out.dir, 'fig-4.jpeg'), width = 1200, height = 600, quality = 100, pointsize = 24)
pdf(paste0(out.dir, 'fig-4.pdf'), width = 16, height = 8)
ggarrange(stress.pca, stress.heatmap[[4]],
          labels = c("A", "B"),
          ncol= 2,
          font.label = list(size = 26))
vp <- viewport(width = 0.5, height = 0.5, x = .41, y = .75)
pushViewport(vp)
print(stress.scree, vp = current.viewport())
popViewport()
dev.off()

#venn
t.stress |>
  mutate(HighCO2 = ifelse(HighCO2>=1, 1, 0), HighLight = ifelse(HighLight>=1, 1, 0), HighSalt = ifelse(HighSalt>=1, 1 , 0), LowN = ifelse(LowN>=1, 1, 0), Heat = ifelse(Heat>=1, 1, 0)) |>
  venn::venn(
    ilabels = 'counts',
    ilcs = 1.8, sncs = 1.8,
    zcolor = cbPalette,
    box = FALSE,
    ggplot = TRUE
  )
ggsave(paste0(supp.outdir, 'venn-stress.jpeg'),
       width = 18, height = 18, dpi = 400)

#high CO2 specific stress proteins: 31
readr::write_tsv(t.stress.sortedco2 |> filter(`High CO2`>=1 & `High Light`<1 & Heat<1 & `High Salt`<1 & `Limited N`<1) |> rownames_to_column(var = 'Geneid'), file.path(analysis.dir, 'tbl.highco2-specific-stress-proteins.tsv'))
#stress proteins that are shared by all five stress conditions - universal stress proteins: 18
readr::write_tsv(t.stress.sortedco2 |> filter(`High CO2`>=1 & `High Light`>=1 & Heat>=1 & `High Salt`>=1 & `Limited N`>=1) |> rownames_to_column(var = 'Geneid'), file.path(analysis.dir, 'tbl.universal-stress-proteins.tsv'))
