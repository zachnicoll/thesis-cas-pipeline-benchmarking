from Bio import Entrez

Entrez.email = 'prospector@qut.com'

n_genomes = 100

print(f"Obtaining {n_genomes} genomes from NCBI...")

rows = (row for row in open('crispr_groundtruth.txt'))
genome_search_terms = (row.split('\t')[4] for row in filter(lambda r: r.find("===") != -1, rows))
ids = [next(genome_search_terms) for _ in range(n_genomes)]

handle = Entrez.efetch(db="nucleotide", id=','.join(ids), rettype="fasta")

output_fasta = open('genomes/genomes.fasta', 'w')
fasta_text = handle.read()
output_fasta.write(fasta_text)

output_fasta.close()
handle.close()

print("Finished obtaining genomes!")
