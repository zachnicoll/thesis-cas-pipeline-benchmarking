# install.packages("rentrez")
# install.packages("tidyverse")

# install.packages("devtools")
# devtools::install_github("ropensci/rentrez")

library(rentrez)
library(tidyverse)
library(stringr)

genome_search_terms = readLines("crispr_groundtruth.txt") %>%
  str_subset(pattern="===") %>%
  str_extract(pattern="\t[ANCPZ_]+[A-Z0-9.]+\t") %>%
  str_replace_all("([\t])", "")


genomes_per_search = 30
max_genomes = 30
found_genomes = 0
fasta = ""

for (i in 1:(length(genome_search_terms)/genomes_per_search)) {
  res <- entrez_search(db = "nuccore", term = paste(genome_search_terms[i*genomes_per_search:i*genomes_per_search + genomes_per_search], collapse=","))
  
  if (length(res$ids) > 0) {
    found_genomes = found_genomes + length(res$ids)
    
    fetched_fasta = entrez_fetch(db = "nuccore", id = res$ids, rettype = "fasta")
    fasta = paste(fasta, fetched_fasta, sep="\n")
  }
  
  if (found_genomes >= max_genomes) {
    break
  }
}

writeLines(fasta, "genomes.fasta")

