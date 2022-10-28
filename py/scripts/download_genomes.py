from Bio import Entrez
import random

Entrez.email = 'prospector@qut.com'

n_genomes = 50

print(f"Obtaining {n_genomes} genomes from NCBI...")
rows = (row for row in open('data_acquisition/crispr_groundtruth.txt'))
genome_search_terms = list((row.split('\t')[4] for row in filter(lambda r: r.find("===") != -1, rows)))

rand_indices = random.sample(range(1, len(genome_search_terms)), n_genomes)
ids = []
for i in rand_indices:
    ids.append(genome_search_terms[i])

handle = Entrez.efetch(db="nucleotide", id=','.join(ids), rettype="fasta")

output_fasta = open('data_acquisition/genomes/genomes.fasta', 'w')
fasta_text = handle.read()
output_fasta.write(fasta_text)

output_fasta.close()
handle.close()

print("Finished obtaining genomes!")
