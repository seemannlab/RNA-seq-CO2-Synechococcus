# Q: Are DE or AA divergent genes over-represented in secretion signals?

library(tidyverse)

# universe of genes to consider
uni <-
  'analysis/I_aa-analysis.tsv' |>
  read_tsv() |>
  pull(Geneid) |>
  unique()

# AA divergent genes
aa.divergent <-
  'analysis/I_aa-analysis.tsv' |>
  read_tsv() |>
  filter(fdr <= .05) |>
  pull(Geneid) |>
  unique() |>
  intersect(uni)

# all DE genes
deg <-
  'analysis/D_stagewise-adjusted-DEGs.tsv' |>
  read_tsv() |>
  filter(is.de) |>
  pull(Geneid) |>
  unique() |>
  intersect(uni)

# genes with secretion signals
secret <-
  'analysis/H_signals.tsv' |>
  read_tsv() |>
  filter(signal != 'No SP') |>
  pull(Geneid) |>
  unique() |>
  intersect(uni)


length(aa.divergent)
# 29
length(deg)
# 2847
length(secret)
# 251

################################################################################

contigency.mat <- function(x, y = secret) {
  # x <- aa.divergent
  # y <- secret
  
  xy <- intersect(x, y) |> length()
  ny <- setdiff(x, y) |> length()
  nx <- setdiff(y, x) |> length()
  nn <- setdiff(uni, union(x, y)) |> length()
  
  # test for greater: 
  # odds for beeing significant (first column) higher for pathway (first row)
  matrix(
    #. Secretion signal,  no signal
    c(xy,                 ny,   # in group of interest X
      nx,                 nn),  # Universe without X
    byrow = TRUE, nrow = 2
  ) |>
    magrittr::set_colnames(c(
      'Secretion signal', 'No signal'
    ))
}

my.p <- function(x) {
  x |>
  fisher.test(alternative = 'greater') |>
  _[['p.value']]
}


################################################################################
de.mat <-
  contigency.mat(deg) |>
  magrittr::set_rownames(c('DE gene', 'not DE'))
de.mat
#              Secretion signal No signal
# DE gene              225      2622
# not DE                26       313

# per chance expected
# 2622 / (313 + 2622) * length(secret)
# 224.2324

my.p(de.mat)
# 0.4918274

################################################################################

aa.mat <-
  contigency.mat(aa.divergent) |>
  magrittr::set_rownames(c('AA divergent', 'not AA divergent'))
aa.mat
#                   Secretion signal No signal
# AA divergent                   10        19
# not AA divergent              241      2916

# per chance expected
# 19 / (2916 + 19) * length(secret)
# 1.624872

my.p(aa.mat)
# 4.031678e-05
