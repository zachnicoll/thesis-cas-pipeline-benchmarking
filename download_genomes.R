# install.packages("rentrez")
# install.packages("tidyverse")

# install.packages("devtools")
# devtools::install_github("ropensci/rentrez")

library(rentrez)
library(tidyverse)
library(stringr)

genome_search_terms = readLines("cas1903.isl.type.txt") %>%
  str_subset(pattern="===") %>%
  str_extract(pattern="\t[ANCPZ_]+[A-Z0-9.]+\t") %>%
  str_replace_all("([\t])", "")

max_genomes = 20
res <- entrez_search(db = "nuccore", term = paste(genome_search_terms[1:max_genomes], collapse=","))
fasta = entrez_fetch(db = "nuccore", id = res$ids, rettype = "fasta")
writeLines(fasta, "genomes.fasta")

