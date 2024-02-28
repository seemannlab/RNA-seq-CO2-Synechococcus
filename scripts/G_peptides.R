# Export peptide sequences of coding genes

library(tidyverse)

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

# exclude all 0s columns
mask <- apply(freqs, 2, max) > 0
freqs2 <- freqs[, mask]
# nicer AA names
colnames(freqs2) <- with(aaMap, set_names(name, let.1))[colnames(freqs2)]

# export freqs
freqs2 |>
  as_tibble(rownames = 'Geneid') |>
  write_tsv('analysis/G_frequencies.tsv')

################################################################################

sessionInfo()
