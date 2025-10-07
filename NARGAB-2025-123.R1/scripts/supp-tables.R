# Collection of longer supplementary tables

################################################################################

library(tidyverse)
library(ggpubr)
library(openxlsx)
library(VennDiagram)
library(GenomicRanges)


library(conflicted)
conflicts_prefer(dplyr::filter)
conflicts_prefer(dplyr::rename)

conflicts_prefer(base::setdiff)

proj.dir <- '/media/dna/projects/rth/co2capture'

################################################################################
# T1: Genomic coordinates

annot <-
  file.path(proj.dir, 'subprojects/RNA-seq_1v15-pct-CO2/RNA-seq-CO2-Synechococcus/data/C_annotation.tsv') |>
  read_tsv()

tbl.refseq <-
  annot |>
  #filter(str_detect(Geneid, '^gene-')) |>
  filter(type != 'rRNA' & type != 'CRS' & type!= 'CRISPR-DR33') |>
  transmute(
    name, type,
    locus_tag, old_locus_tag,
    product,
    chromosome = seqnames,
    start, end, strand
  ) |>
  arrange(chromosome, start)

################################################################################
# T2: RNA-seq mapping stats

tbl.rnaseq <-
  file.path(proj.dir, 'subprojects/RNA-seq_1v15-pct-CO2/RNA-seq-CO2-Synechococcus/data/A_general-stat.tsv') |>
  read_tsv()


################################################################################
# T3: TPM

tbl.tpm <-
  file.path(proj.dir, 'continuation-ses/subprojects/RNA-seq_1v15-pct-CO2/RNA-seq-CO2-Synechococcus/analysis/D_tpm-expression.tsv') |>
  read_tsv()


################################################################################
# T4: average z-scaled gene expression

tbl.zscale <-
  file.path(proj.dir, 'continuation-ses/subprojects/RNA-seq_1v15-pct-CO2/RNA-seq-CO2-Synechococcus/analysis/D_z-expr-mean-cond.tsv') |>
  read_tsv() |>
  rename('0.04% CO2' = '0.04',
         '4% CO2' = '4',
         '8% CO2' = '8',
         '30% CO2' = '30')


################################################################################
# T5: Number of differential expressed genes (DEGs; DESeq2 FDR ≤ 0.05) for different log2 fold changes

deg.stats <- file.path(proj.dir, 'continuation-ses/subprojects/RNA-seq_1v15-pct-CO2/RNA-seq-CO2-Synechococcus/analysis/D_overview-by-logFC.tsv') |>
  read_tsv()
tbl.deg.stats <-
  deg.stats |>
  dplyr::rename(comparison = test)


################################################################################
# T6:Differential expression analysis and annotations

deg <-
  file.path(proj.dir, 'continuation-ses/subprojects/RNA-seq_1v15-pct-CO2/RNA-seq-CO2-Synechococcus/analysis/D_DEGs.tsv') |>
  read_tsv()

tbl.deg <-
  deg |>
  arrange(desc(baseMean), Geneid, is.de, log2FoldChange) |>
  dplyr::rename(
    comparison = test,
    'average expression across all conditions' = baseMean,
    'adjusted p-value (FDR)' = padj,
    'diff. expressed at FDR<=0.05' = is.de
  ) |>
  select(-pvalue, - seqnames, - start, - end, - strand) |> 
  select(Geneid, type, name, locus_tag, old_locus_tag, product, everything()) |>
  mutate_if(is.character, replace_na, '')


################################################################################
# T7: Venn diagram of DEGs (Figure 1C)

deg |> mutate(log2FoldChange = round(log2FoldChange, 2)) |> filter(is.de & abs(log2FoldChange) >= 1) |> group_by(test) |> group_split() |> as.list() -> foo
set_names(
  foo |> map(pull, Geneid),
  foo |> map(pull, test) |> map(dplyr::first) |> unlist()
) |>
  _[c(
    '0.04->4+8% CO2', '4+8->30% CO2',
    '0.04->30% CO2', 'Other->30% CO2'
  )] |>
  set_names(
    '0.04% vs 4+8%',
    '4+8% vs 30%',
    '0.04% vs 30%',
    'Other vs 30%'
  ) -> de.list
deg.venn <- VennDiagram::get.venn.partitions(de.list)
tbl.venn <-
  deg.venn |>
  select(`0.04% vs 4+8%`, `4+8% vs 30%`, `0.04% vs 30%`, `Other vs 30%`, ..values.., ..count..) |>
  rename("Geneid" = ..values.., "Count" = ..count..) |>
  unnest(Geneid)


################################################################################
# T8: gene to KEGG pathway matches

tbl.kegg <-
  file.path(proj.dir, 'subprojects/RNA-seq_1v15-pct-CO2/RNA-seq-CO2-Synechococcus/analysis/H_gene2pathway.tsv') |>
  read_tsv() |>
  select(
    pathway = Pathway,
    title = Title,
    gene= old_locus_tag,
    name
  )


################################################################################
# T9: GSEA in KEGG pathways using log2FC (Figure 3A, S5, and S6)

kegg.1.1 <- tbl.kegg %>% filter(pathway == "syp00195") %>% rename('old_locus_tag' = 'gene') %>% left_join(tbl.deg, by = 'old_locus_tag') %>% arrange(desc(log2FoldChange)) %>% mutate(comparison = "Other->30% CO2")
kegg.1.2 <- tbl.kegg %>% filter(pathway == "syp00190") %>% rename('old_locus_tag' = 'gene') %>% left_join(tbl.deg, by = 'old_locus_tag') %>% arrange(desc(log2FoldChange)) %>% mutate(comparison = "Other->30% CO2")
kegg.1.3 <- tbl.kegg %>% filter(pathway == "syp00196") %>% rename('old_locus_tag' = 'gene') %>% left_join(tbl.deg, by = 'old_locus_tag') %>% arrange(desc(log2FoldChange)) %>% mutate(comparison = "Other->30% CO2")
kegg.1.4 <- tbl.kegg %>% filter(pathway == "syp00860") %>% rename('old_locus_tag' = 'gene') %>% left_join(tbl.deg, by = 'old_locus_tag') %>% arrange(desc(log2FoldChange)) %>% mutate(comparison = "Other->30% CO2")

kegg.2.1 <- tbl.kegg %>% filter(pathway == "syp03010") %>% rename('old_locus_tag' = 'gene') %>% left_join(tbl.deg, by = 'old_locus_tag') %>% arrange(desc(log2FoldChange)) %>% mutate(comparison = "0.04->30% CO2")
kegg.2.2 <- tbl.kegg %>% filter(pathway == "syp00920") %>% rename('old_locus_tag' = 'gene') %>% left_join(tbl.deg, by = 'old_locus_tag') %>% arrange(desc(log2FoldChange)) %>% mutate(comparison = "0.04->30% CO2")
kegg.2.3 <- tbl.kegg %>% filter(pathway == "syp00190") %>% rename('old_locus_tag' = 'gene') %>% left_join(tbl.deg, by = 'old_locus_tag') %>% arrange(desc(log2FoldChange)) %>% mutate(comparison = "0.04->30% CO2")
kegg.2.4 <- tbl.kegg %>% filter(pathway == "syp00790") %>% rename('old_locus_tag' = 'gene') %>% left_join(tbl.deg, by = 'old_locus_tag') %>% arrange(desc(log2FoldChange)) %>% mutate(comparison = "0.04->30% CO2")
kegg.2.5 <- tbl.kegg %>% filter(pathway == "syp01212") %>% rename('old_locus_tag' = 'gene') %>% left_join(tbl.deg, by = 'old_locus_tag') %>% arrange(desc(log2FoldChange)) %>% mutate(comparison = "0.04->30% CO2")
kegg.2.6 <- tbl.kegg %>% filter(pathway == "syp00400") %>% rename('old_locus_tag' = 'gene') %>% left_join(tbl.deg, by = 'old_locus_tag') %>% arrange(desc(log2FoldChange)) %>% mutate(comparison = "0.04->30% CO2")
kegg.2.7 <- tbl.kegg %>% filter(pathway == "syp00195") %>% rename('old_locus_tag' = 'gene') %>% left_join(tbl.deg, by = 'old_locus_tag') %>% arrange(desc(log2FoldChange)) %>% mutate(comparison = "0.04->30% CO2")

kegg.3.1 <- tbl.kegg %>% filter(pathway == "syp03010") %>% rename('old_locus_tag' = 'gene') %>% left_join(tbl.deg, by = 'old_locus_tag') %>% arrange(desc(log2FoldChange)) %>% mutate(comparison = "4+8->30% CO2")
kegg.3.2 <- tbl.kegg %>% filter(pathway == "syp00195") %>% rename('old_locus_tag' = 'gene') %>% left_join(tbl.deg, by = 'old_locus_tag') %>% arrange(desc(log2FoldChange)) %>% mutate(comparison = "4+8->30% CO2")
kegg.3.3 <- tbl.kegg %>% filter(pathway == "syp00196") %>% rename('old_locus_tag' = 'gene') %>% left_join(tbl.deg, by = 'old_locus_tag') %>% arrange(desc(log2FoldChange)) %>% mutate(comparison = "4+8->30% CO2")
kegg.3.4 <- tbl.kegg %>% filter(pathway == "syp00860") %>% rename('old_locus_tag' = 'gene') %>% left_join(tbl.deg, by = 'old_locus_tag') %>% arrange(desc(log2FoldChange)) %>% mutate(comparison = "4+8->30% CO2")
kegg.3.5 <- tbl.kegg %>% filter(pathway == "syp00970") %>% rename('old_locus_tag' = 'gene') %>% left_join(tbl.deg, by = 'old_locus_tag') %>% arrange(desc(log2FoldChange)) %>% mutate(comparison = "4+8->30% CO2")
kegg.3.6 <- tbl.kegg %>% filter(pathway == "syp00910") %>% rename('old_locus_tag' = 'gene') %>% left_join(tbl.deg, by = 'old_locus_tag') %>% arrange(desc(log2FoldChange)) %>% mutate(comparison = "4+8->30% CO2")
kegg.3.7 <- tbl.kegg %>% filter(pathway == "syp00190") %>% rename('old_locus_tag' = 'gene') %>% left_join(tbl.deg, by = 'old_locus_tag') %>% arrange(desc(log2FoldChange)) %>% mutate(comparison = "4+8->30% CO2")
kegg.3.8 <- tbl.kegg %>% filter(pathway == "syp00220") %>% rename('old_locus_tag' = 'gene') %>% left_join(tbl.deg, by = 'old_locus_tag') %>% arrange(desc(log2FoldChange)) %>% mutate(comparison = "4+8->30% CO2")

tbl.kegg.lfc <- rbind(kegg.1.1, kegg.1.2, kegg.1.3, kegg.1.4,
                      kegg.2.1, kegg.2.2, kegg.2.3, kegg.2.4, kegg.2.5, kegg.2.6, kegg.2.7,
                      kegg.3.1, kegg.3.2, kegg.3.3, kegg.3.4, kegg.3.5, kegg.3.6, kegg.3.7, kegg.3.8) |>
  select(comparison, pathway, title, Geneid, product, log2FoldChange, `adjusted p-value (FDR)`)


################################################################################
# T10 and T11: STRING network clustering and enrichment

mcl <- 
  file.path(proj.dir, 'continuation-ses/subprojects/RNA-seq_1v15-pct-CO2/RNA-seq-CO2-Synechococcus/analysis/K_string-mcl-clustering.csv') |>
  read_delim(delim = ",") |>
  filter(! is.na(`__mclCluster`)) |>
  transmute(
    cluster = paste0('cluster_', `__mclCluster`),
    string.name = `display name`,
    locus = `query term`,
    product = `stringdb::description`
  ) |>
  arrange(cluster, string.name)

# filter to ≥3 genes
string.mcl <-
  mcl |>
  dplyr::count(cluster) |>
  filter(n >= 3) |>
  arrange(desc(n)) |>
  select(cluster) |>
  left_join(mcl, 'cluster')

string.enrich <- 
  file.path(proj.dir, 'continuation-ses/subprojects/RNA-seq_1v15-pct-CO2/RNA-seq-CO2-Synechococcus/analysis/K_string-enrichment.csv') |>
  read_delim(delim = ",") |>
  transmute(
    cluster = paste0('cluster_', cluster),
    database = category,
    geneset = `term name`,
    `geneset description` = description,
    `# genes`,
    `# background genes`,
    `p-value`,
    FDR = `FDR value`,
    `gene list` = genes,
  )


################################################################################
# T12: Stress proteins

tbl.stress <- file.path(proj.dir, 'continuation-ses/subprojects/RNA-seq_1v15-pct-CO2/RNA-seq-CO2-Synechococcus/stress/tbl.stress.tsv') |> read_tsv()


################################################################################
# T13: orphan RNAs

tbl.orphanrna <- 
  file.path(proj.dir, 'continuation-ses/subprojects/RNA-seq_1v15-pct-CO2/RNA-seq-CO2-Synechococcus/orphan-rna/featurecounts-readpair/intergenic-offset100-nostrand-featurecounts-readpairs-l100-c5kpercond-annotated.txt') |>
  read_tsv(col_names = F) |>
  select(
    `Adjacent genes` = X1,
    'Chromosome' = X2,
    'Start' = X3,
    'End' = X4,
    'Strand' = X5,
    'Length' = X6,
    'S1' = X7,
    'S2' = X8,
    'S3' = X9,
    'S4' = X10,
    'S5' = X11,
    'S6' = X12,
    'S7' = X13,
    'S8' = X14,
    'S9' = X15,
    'S10' = X16,
    'S11' = X17,
    'S12' = X18,
    'S13' = X19,
    'S14' = X20,
    'S15' = X21,
    'S16' = X22,
    'Annotation' = X24
  )

#get TPM
tbl.raw.noribo <-
  file.path(proj.dir, 'continuation-ses/subprojects/RNA-seq_1v15-pct-CO2/RNA-seq-CO2-Synechococcus/analysis/D_raw-expression-noribo.tsv') |>
  read_tsv()

tbl.raw.all <- 
  tbl.orphanrna |>
  select(
    'Geneid' = `Adjacent genes`,
    S1, S2, S3, S4, S5, S6, S7, S8, S9, S10, S11, S12, S13, S14, S15, S16,
    'length' = Length
  ) |>
  bind_rows(tbl.raw.noribo)

library(DGEobj.utils)
#see https://haroldpimentel.wordpress.com/2014/05/08/what-the-fpkm-a-review-rna-seq-expression-units/
tpm <- convertCounts(tbl.raw.all[,-ncol(tbl.raw.all)] |> column_to_rownames('Geneid') |> as.matrix(), unit = "TPM", geneLength = tbl.raw.all[,ncol(tbl.raw.all)]$length)
tpm <- tpm %>% as.data.frame() %>%
  rownames_to_column("rowname") %>%
  rowwise() %>% mutate(min_tpm = min(c_across(-rowname)), max_tpm = max(c_across(-rowname)))
write.xlsx(
  tpm,
  file.path(proj.dir, 'continuation-ses/subprojects/RNA-seq_1v15-pct-CO2/RNA-seq-CO2-Synechococcus/analysis/tbl.tpm-expression-orphan.xlsx'),
  colWidths = 'auto'
)

tbl.orphanrna <- 
  tbl.orphanrna |>
  select(`Adjacent genes`, Chromosome, Start, End, Strand, Length, Annotation) |>
  left_join(tpm |> rename(`Adjacent genes` = 'rowname'), 'Adjacent genes') |>
  rename('min TPM' = 'min_tpm', 'max TPM' = 'max_tpm')

#log2FC_0.04->30 = log2 (mean(30%)/mean(0.04%))
tbl.orphanrna <- tbl.orphanrna |> mutate(`log2FC 0.04->30%` = log2((S13+S14+S15+S16)/(S1+S2+S3+S4)))
#log2FC_4->30 = log2 (mean(30%)/mean(4%+8%))
tbl.orphanrna <- tbl.orphanrna |> mutate(`log2FC 4+8->30%` = log2(((S13+S14+S15+S16)/4)/((S5+S6+S7+S8+S9+S10+S11+S12)/8)))


################################################################################
# T14: Cisreg RNAs (CRSs and riboswitches) and pseudogenes

gff <- file.path(proj.dir, 'continuation-ses/subprojects/RNA-seq_1v15-pct-CO2/RNA-seq-CO2-Synechococcus/featurecounts-no-genes/refseq-no-genes.gff') |>
  rtracklayer::import.gff3()

annot.gff <- (gff |> mcols())[,c("source", "type", "ID", "Name", "Note")] |> as.data.frame() |> rename("Geneid" = ID)

tbl.cisreg.pseudo <-
  file.path(proj.dir, 'continuation-ses/subprojects/RNA-seq_1v15-pct-CO2/RNA-seq-CO2-Synechococcus/featurecounts-no-genes/refseq-no-genes-featurecounts-readpairs.txt') |>
  read_tsv(col_names = T) |>
  select(
    Geneid,
    'Chromosome' = Chr,
    Start, End, Strand, Length,
    'S1' = 'CO2-0-04percent_Lane1-S1-featurecounts.txt',
    'S2' = 'CO2-0-04percent_Lane2-S2-featurecounts.txt',
    'S3' = 'CO2-0-04percent_Lane2-S3-featurecounts.txt',
    'S4' = 'CO2-0-04percent_Lane2-S4-featurecounts.txt',
    'S5' = 'CO2-4percent_Lane2-S5-featurecounts.txt',
    'S6' = 'CO2-4percent_Lane1-S6-featurecounts.txt',
    'S7' = 'CO2-4percent_Lane1-S7-featurecounts.txt',
    'S8' = 'CO2-4percent_Lane2-S8-featurecounts.txt',
    'S9' = 'CO2-8percent_Lane2-S9-featurecounts.txt',
    'S10' = 'CO2-8percent_Lane2-S10-featurecounts.txt',
    'S11' = 'CO2-8percent_Lane2-S11-featurecounts.txt',
    'S12' = 'CO2-8percent_Lane1-S12-featurecounts.txt',
    'S13' = 'CO2-30percent_Lane2-S13-featurecounts.txt',
    'S14' = 'CO2-30percent_Lane1-S14-featurecounts.txt',
    'S15' = 'CO2-30percent_Lane2-S15-featurecounts.txt',
    'S16' = 'CO2-30percent_Lane1-S16-featurecounts.txt'
  )

tbl.raw.all <- 
  tbl.cisreg.pseudo |>
  select(
    Geneid, S1, S2, S3, S4, S5, S6, S7, S8, S9, S10, S11, S12, S13, S14, S15, S16,
    'length' = Length
  ) |>
  bind_rows(tbl.raw.noribo)

library(DGEobj.utils)
#see https://haroldpimentel.wordpress.com/2014/05/08/what-the-fpkm-a-review-rna-seq-expression-units/
tpm <- convertCounts(tbl.raw.all[,-ncol(tbl.raw.all)] |> column_to_rownames('Geneid') |> as.matrix(), unit = "TPM", geneLength = tbl.raw.all[,ncol(tbl.raw.all)]$length)
tpm <- tpm %>% as.data.frame() %>%
  rownames_to_column("rowname") %>%
  rowwise() %>% mutate(min_tpm = min(c_across(-rowname)), max_tpm = max(c_across(-rowname)))
write.xlsx(
  tpm,
  file.path(proj.dir, 'continuation-ses/subprojects/RNA-seq_1v15-pct-CO2/RNA-seq-CO2-Synechococcus/analysis/tbl.tpm-expression-cisreg-pseudo.xlsx'),
  colWidths = 'auto'
)

tbl.cisreg.pseudo <- 
  annot.gff |>
  left_join(tpm |> rename('Geneid' = 'rowname'), 'Geneid') |>
  left_join(tbl.cisreg.pseudo |> select(Geneid, Chromosome, Start, End, Strand, Length), 'Geneid') |>
  select(Geneid, Source = source, Type = type, Name, Note, Chromosome, Start, End, Strand, Length, S1, S2, S3, S4, S5, S6, S7, S8, S9, S10, S11, S12, S13, S14, S15, S16, 'min TPM' = 'min_tpm', 'max TPM' = 'max_tpm')

#log2FC_0.04->30 = log2 (mean(30%)/mean(0.04%))
tbl.cisreg.pseudo <- tbl.cisreg.pseudo |> mutate(`log2FC 0.04->30%` = ifelse(!is.na(log2((S13+S14+S15+S16+1)/(S1+S2+S3+S4+1))), log2((S13+S14+S15+S16+1)/(S1+S2+S3+S4+1)), 0))
#log2FC_4->30 = log2 (mean(30%)/mean(4%+8%))
tbl.cisreg.pseudo <- tbl.cisreg.pseudo |> mutate(`log2FC 4+8->30%` = ifelse(!is.na(log2(((S13+S14+S15+S16)/4)/((S5+S6+S7+S8+S9+S10+S11+S12)/8))), log2(((S13+S14+S15+S16)/4)/((S5+S6+S7+S8+S9+S10+S11+S12)/8)), 0))


################################################################################
# T15: AA-composition

tbl.aa <-
  file.path(proj.dir, 'subprojects/RNA-seq_1v15-pct-CO2/RNA-seq-CO2-Synechococcus/analysis/I_aa-analysis.tsv') |>
  read_tsv() |>
  filter(fdr <= 0.05) |>
  left_join(tbl.deg, 'Geneid') |>
  filter(comparison == "Other->30% CO2") |>
  select(
    Geneid, name, old_locus_tag, product, AA,
    Observed = value,
    Length = length.x,
    Expected = expected,
    `FDR AA` = fdr,
    `log2FoldChange Other vs 30% CO2` = log2FoldChange,
    `FDR Other vs 30% CO2` = `adjusted p-value (FDR)`
  )

write.xlsx(
  tbl.aa,
  file.path(proj.dir, 'continuation-ses/subprojects/RNA-seq_1v15-pct-CO2/RNA-seq-CO2-Synechococcus/analysis/tbl.aa-composition.xlsx'),
  colWidths = 'auto'
)

tbl.aa %>% filter(`FDR Other vs 30% CO2`<=0.05) %>% distinct(Geneid) %>% dplyr::count()


################################################################################
# T16: SignalP enrichment and results

tbl.signal <-
  file.path(proj.dir, 'subprojects/RNA-seq_1v15-pct-CO2/RNA-seq-CO2-Synechococcus/analysis/H_signals.tsv') |>
  read_tsv() |>
  select(
    locus_tag, old_locus_tag,
    name,
    signal,
    chromosome = seqnames,
    start, end, strand,
  ) |>
  arrange(chromosome, start)


# additional data retrieval about SignalP
tbl.signal.genes <-
  file.path(proj.dir, 'continuation-ses/subprojects/RNA-seq_1v15-pct-CO2/RNA-seq-CO2-Synechococcus/analysis/H_signal-enrichment-genes.tsv') |>
  read_tsv() |>
  rename('old_locus_tag' = 'genes') |>
  left_join(tbl.deg, 'old_locus_tag')

#24 out of 251 genes with predicted secretion signal are differential expressed in Others vs 30% CO2 (Fisher Exact Test p=0.06414)
tbl.signal %>% filter(signal != "No SP") %>% left_join(tbl.deg, 'locus_tag') %>% filter(comparison == "Other->30% CO2" & (`adjusted p-value (FDR)` <= 0.05 & abs(log2FoldChange)>=1)) %>% dplyr::count()
deg.sp <- matrix(c(24, 251-24, 215, 3186-215), nrow = 2)
fisher.test(deg.sp, alternative = "greater")

#87 of 251 genes with predicted secretion signal are differential expressed in at least one pairwise CO2 comparison (Fisher Exact Test p=0.3536)
tbl.signal %>% filter(signal != "No SP") %>% left_join(tbl.deg, 'locus_tag') %>% filter(`adjusted p-value (FDR)` <= 0.05 & abs(log2FoldChange)>=1) %>% distinct(locus_tag) %>% dplyr::count()
deg <- matrix(c(87, 251-87, 1061, 3186-1061), nrow = 2)
fisher.test(deg, alternative = "greater")


################################################################################
export <- list(
  'S1 RefSeq genes + Rfam sRNAs' = tbl.refseq, #add tbl.rfam
  'S2 RNA-seq statistics' = tbl.rnaseq,
  'S3 TPM-normalized counts' = tbl.tpm,
  'S4 Average z-scaled expression' = tbl.zscale,
  'S5 DEGs statistics' = tbl.deg.stats,
  'S6 DEGs' = tbl.deg,
  'S7 Venn DEGs' = tbl.venn,
  'S8 KEGG pathways' = tbl.kegg,
  'S9 KEGG GSEA of DEGs' = tbl.kegg.lfc,
  'S10 Network clustering' = string.mcl,
  'S11 Network enrichment' = string.enrich,
  'S12 Stress proteins' = tbl.stress,
  'S13 Orphan RNAs' = tbl.orphanrna,
  'S14 Pseudogene+Riboswitch+CRS' = tbl.cisreg.pseudo
  #'S15 Amino acid composition' = tbl.aa,
  #'S16 SignalP' = tbl.signal
)


################################################################################
################################################################################
# Export excel style table

hs <- createStyle(
  textDecoration = "BOLD", fontColour = "#FFFFFF", fontSize = 15,
  fontName = "Arial Narrow", fgFill = "#4F80BD"
)

bs <- createStyle(
  fontColour = "#000000", fontSize = 16,
  fontName = "Arial Narrow", fgFill = "#FFFFFF"
)

write.xlsx(
  export,
  file.path(proj.dir, 'continuation-ses/subprojects/RNA-seq_1v15-pct-CO2/Manuscript-Overleaf/RNA-CO2-supp-tables.xlsx'),
  colWidths = 'auto',
  headerStyle = hs,
  tableStyle = bs
)

################################################################################
################################################################################

