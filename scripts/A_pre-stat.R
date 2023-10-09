#!/usr/bin/env Rscript

# Create a summary table of the mapping/filtering/etc statistics

library(tidyverse)

################################################################################

# Helper to quickly load a summary of the multiQC json
load.json <- function(i) {
  dat <- jsonlite::read_json(i)
  
  dat$report_general_stats_data %>%
    map(function(xs) {
      map2(xs, names(xs),
           function(x, y) {
             x$sample <- y
             keep(x, ~ length(.x) == 1) %>%
               as_tibble()
           }) %>%
        bind_rows()
    }) %>%
    bind_rows()
}
# multiqc might add suffixes -> remove and collect non-NA entries
combine.sample <- function(xs) {
  xs %>%
    mutate_at('sample', str_remove, '-fastp$') %>%
    mutate_at('sample', str_remove, '_R2$') %>%
    group_by(sample) %>%
    summarize_all(discard, is.na) %>%
    ungroup
}

################################################################################


c(
  './raw-data/schlange-data1-multiqc.json',
  './raw-data/schlange-data2-multiqc.json',
  './raw-data/public-datasets/multiqc-schlange-illumina.json',
  './raw-data/public-datasets/multiqc-schlange-solid.json'
) %>%
  map(function(i) {
    i %>%
      load.json() %>%
      select(
        sample,
        before_filtering_total_reads,
        after_filtering_total_reads,
        non_rRNA
      ) -> xs
    # make numbers to pairs for paired-end
    if(any(str_detect(xs$sample, '_R2$'))) {
      xs %>%
        mutate_if(is.numeric, ~ .x / 2) -> xs
    }
    xs %>%
      combine.sample() %>%
      separate(sample, c('batch', 'sample', 'condition'), sep = '_')
    
    
  }) %>%
  bind_rows() -> schlange

################################################################################

c(
  './raw-data/browser-data1-multiqc.json',
  './raw-data/browser-data2-multiqc.json',
  './raw-data/public-datasets/multiqc-public-browser.json'
) %>%
  map(function(i) {
    i %>%
      load.json() -> foo
    
    if('unpaired_aligned_none' %in% colnames(foo)) {
      foo$paired_aligned_none <- foo$unpaired_aligned_none
    }
    
    foo %>%
      transmute(
        sample,
        aligned = total_reads - paired_aligned_none,
        filtered.aln = Total / 2,
        Assigned = Assigned / 2,
      ) %>%
      combine.sample() %>%
      separate(sample, c('condition', 'sample', 'dataset'), sep = '_')
  }) %>%
  bind_rows() -> browser

################################################################################
# Combine 

dat <-
  full_join(
    schlange,
    browser %>% mutate_at('sample', str_remove, '^Lane[0-9]-'),
    'sample'
  )

################################################################################
# Nice table for the CO2 libs of this project

res.co2 <-
  dat |>
  filter(! str_detect(sample, '^SRR[0-9]+$')) |>
  transmute(
    'CO2%' = condition.y %>%
      str_remove('^CO2-') %>%
      str_remove('percent') %>%
      str_replace('-', '.') %>%
      as.numeric,
    sample,
    batch,
    Reads = prettyNum(before_filtering_total_reads, big.mark = ','),
    'Cleaned reads' = sprintf(
      '%s (%.1f%%)',
      prettyNum(after_filtering_total_reads, big.mark = ','), 
      after_filtering_total_reads / before_filtering_total_reads * 100
    ),
    'Non rRNA' = sprintf(
      '%s (%.1f%%)',
      prettyNum(non_rRNA, big.mark = ','), 
      non_rRNA / after_filtering_total_reads * 100
    ),
    'Mapped' = sprintf(
      '%s (%.1f%%)',
      prettyNum(aligned, big.mark = ','), 
      aligned / non_rRNA * 100
    ),
    'Filtered alignment' = sprintf(
      '%s (%.1f%%)',
      prettyNum(filtered.aln, big.mark = ','), 
      filtered.aln / aligned * 100
    ),
    'Assigned' = sprintf(
      '%s (%.1f%%)',
      prettyNum(Assigned, big.mark = ','), 
      Assigned / filtered.aln * 100
    )
  ) %>%
  arrange(`CO2%`, sample)

write_tsv(res.co2, 'data/A_general-stat.tsv')

################################################################################
# Nice table for the public datasets

'./raw-data/public-datasets/multiqc-schlange-solid.json' |>
  load.json() |>
  pull(sample) |>
  str_extract('SRR[0-9]*') -> solid.libs


res.public <-
  dat |>
  filter(str_detect(sample, '^SRR[0-9]+$')) |>
  transmute(
    Platform = ifelse(sample %in% solid.libs, 'SOLiD', 'Illumina'),
    Dataset = dataset,
    Condition = condition.y |>
      str_remove('^replicate.-'),
    SRA = sample,
    Reads = prettyNum(before_filtering_total_reads, big.mark = ','),
    'Cleaned reads' = sprintf(
      '%s (%.1f%%)',
      prettyNum(after_filtering_total_reads, big.mark = ','), 
      after_filtering_total_reads / before_filtering_total_reads * 100
    ),
    'Non rRNA' = sprintf(
      '%s (%.1f%%)',
      prettyNum(non_rRNA, big.mark = ','), 
      non_rRNA / after_filtering_total_reads * 100
    ),
    'Mapped' = sprintf(
      '%s (%.1f%%)',
      prettyNum(aligned, big.mark = ','), 
      aligned / non_rRNA * 100
    ),
    'Filtered alignment' = sprintf(
      '%s (%.1f%%)',
      prettyNum(filtered.aln, big.mark = ','), 
      filtered.aln / aligned * 100
    ),
    'Assigned' = sprintf(
      '%s (%.1f%%)',
      prettyNum(Assigned, big.mark = ','), 
      Assigned / filtered.aln * 100
    )
  ) |>
  arrange(Platform, Dataset, Condition, SRA)

write_tsv(res.public, 'data/A_public-stat.tsv')
