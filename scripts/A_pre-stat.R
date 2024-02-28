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
# QC processing stats

schlange <-
  'raw-data/schlange-data-multiqc.json' |>
  load.json() |>
  select(
    sample,
    before_filtering_total_reads,
    after_filtering_total_reads,
    non_rRNA
  ) |>
  combine.sample() |>
  separate(sample, c('batch', 'sample', 'condition'), sep = '_')

################################################################################
# Mapping stats

browser <-
  'raw-data/browser-data-multiqc.json' |>
  load.json() |>
  transmute(
    sample,
    aligned = total_reads - paired_aligned_none,
    filtered.aln = Total / 2,
    Assigned = Assigned / 2,
  ) |>
  combine.sample() |>
  separate(sample, c('condition', 'sample', 'dataset'), sep = '_')

################################################################################
# Combine 

dat <-
  full_join(
    schlange,
    browser |> mutate_at('sample', str_remove, '^Lane[0-9]-'),
    'sample'
  )

################################################################################
# Nice summary table

res.co2 <-
  dat |>
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
sessionInfo()
