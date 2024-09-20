library(tidyverse)

################################################################################
# 'normal' annotation

tbl.annot <-
  'data/C_annotation.tsv' |>
  read_tsv()

################################################################################
# averaged zscaled expression

tbl.z <-
  'analysis/E_zvst.tsv' |>
  read_tsv() |>
  select(Geneid, CO2, avg.z) |>
  mutate(CO2 = sprintf('Z-Scaled expression %s%%', CO2)) |>
  spread(CO2, avg.z) |>
  select(
    Geneid,
    `Z-Scaled expression 0.04%`,
    `Z-Scaled expression 4%`,
    `Z-Scaled expression 8%`,
    `Z-Scaled expression 30%`
  )

################################################################################
# ugly 16x2 deg part

tbl.deg <-
  'analysis/D_stagewise-adjusted-DEGs.tsv' |>
  read_tsv() |>
  select(Geneid, test, log2FoldChange, padj.stageR) |>
  pivot_longer(c(log2FoldChange, padj.stageR)) |>
  mutate(
    x = paste(
      test |> str_remove(' CO2$'),
      ifelse(name == 'log2FoldChange', 'log2FC', 'FDR')
    )
  ) |>
  select(Geneid, name = x, value) |>
  pivot_wider()

# indication in at least one pairwise comparison DE
tbl.dextra <-
  'analysis/D_stagewise-adjusted-DEGs.tsv' |>
  read_tsv() |>
  filter(is.de) |>
  transmute(
    Geneid,
    'DE in ≥1 pairwsie comparisons' = 'yes'
  ) |>
  unique()

################################################################################
# Something KEGG
tbl.brite <-
  'analysis/E_brite-hierarchy.tsv' |>
  read_tsv() |>
  filter(level1 == 'KEGG Orthology (KO)') |>
  select(- c(level1, level5, level6)) |>
  mutate(row = 1:n()) |>
  pivot_longer(- c(old_locus_tag, row)) |>
  group_by(old_locus_tag, row) |>
  summarize(brite = value |>
              discard(is.na) |>
              unlist() |>
              str_c(collapse = ' > ') |>
              str_remove(' > $')) |>
  group_by(old_locus_tag) |>
  mutate(v = 1:n()) |>
  ungroup() |>
  transmute(
    old_locus_tag,
    name = paste('KEGG orthology', v),
    value = brite
  ) |>
  pivot_wider()

################################################################################
# KEGG pathways

tbl.path <-
  'analysis/H_gene2pathway.tsv' |>
  read_tsv() |>
  transmute(
    Geneid,
    value = sprintf('%s (%s)', Title, Pathway)
  ) |>
  group_by(Geneid) |>
  mutate(name = paste('KEGG pathway', 1:n())) |>
  pivot_wider()


################################################################################
# ugly go

tbl.monster <-
  tbl.annot |>
  left_join(tbl.z, 'Geneid') |>
  left_join(tbl.dextra, 'Geneid') |>
  mutate_at('DE in ≥1 pairwsie comparisons', replace_na, 'no') |>
  left_join(tbl.deg, 'Geneid') |>
  left_join(tbl.brite, 'old_locus_tag') |>
  left_join(tbl.path, 'Geneid') 

write_tsv(tbl.monster, 'nuf-extra-monster.tsv')

################################################################################
# counts separately

tbl.annot |>
  inner_join(
    'analysis/D_normalized-counts.tsv' |>
      read_tsv() |>
    select(
      `Geneid`,
      `CO2-0-04percent_Lane1-S1_CO2-adaptation`,
      `CO2-0-04percent_Lane2-S2_CO2-adaptation`, 
      `CO2-0-04percent_Lane2-S3_CO2-adaptation`,
      `CO2-0-04percent_Lane2-S4_CO2-adaptation`, 
      `CO2-4percent_Lane2-S5_CO2-adaptation`,
      `CO2-4percent_Lane1-S6_CO2-adaptation`,
      `CO2-4percent_Lane1-S7_CO2-adaptation`, 
      `CO2-4percent_Lane2-S8_CO2-adaptation`, 
      `CO2-8percent_Lane2-S9_CO2-adaptation`,
      `CO2-8percent_Lane2-S10_CO2-adaptation`, 
      `CO2-8percent_Lane2-S11_CO2-adaptation`,
      `CO2-8percent_Lane1-S12_CO2-adaptation`,
      `CO2-30percent_Lane2-S13_CO2-adaptation`,
      `CO2-30percent_Lane1-S14_CO2-adaptation`,
      `CO2-30percent_Lane2-S15_CO2-adaptation`, 
      `CO2-30percent_Lane1-S16_CO2-adaptation`
    ),
    'Geneid'
  ) |>
  write_tsv('nuf-counts.tsv')

################################################################################
