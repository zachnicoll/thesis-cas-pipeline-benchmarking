============================================================
===	Supplementary Dataset 1
===	Supplementary_Dataset_1.loci.tgz
===	Annotated CRISPR-Cas loci
============================================================
The gzipped tar archive contains four plain ASCII tab-delimited text files:

cas1903.isl.type.txt	- complete, fully annotated single-system loci

cas1903.isl.mult.txt	- loci, contaning several merged or adjacent CRISPR-Cas systems, where each system is complete and fully annotated

cas1903.isl.part.txt	- loci, contaning incomplete systems

cas1903.isl.typx.txt	- loci, contaning unclassified systems

Prok1903.gen.txt	- complete list of genomes in the dataset

Each of the loci files contains a list of loci; each locus consists of a header line and a number of lines, describing genes.

............................................................
...	Loci header line
............................................................
1	- "===", header line tag

2	- locus ID

3	- number of genes in locus

4	- NCBI genome assembly ID

5	- GenBank ID of the genome partition/contig sequences

6-10	- blank

11	- system (sub)type

12	- genome name (non-unique; different assemblies can have the same name)

............................................................
...	Loci gene description
............................................................
1	- gene ID (NCBI locus tag or custom unique ID)

2	- coordinates in the nucleotide sequence

3	- coding direction

4	- NCBI genome assembly ID

5	- GenBank ID of the genome partition/contig sequences

6	- ordinal number of gene in the genome partition/contig sequences

7	- local protein sequence ID (derived from GenBank accession-dot-version ID, e.g. WP_096410491_1 is derived from WP_096410491.1)

8	- "+" marks CRISPR-Cas genes

9	- profile(s) detected in the protein sequence (comma-separated list)

10	- sequence family annotation (multi-domain proteins are annotated as a comma-separated list of IDs)

11	- system (sub)type

12	- genome name (non-unique; different assemblies can have the same name)

............................................................
...	Genome file
............................................................
1	- NCBI genome assembly ID

2	- genome assembly weight

3	- genome name (non-unique; different assemblies can have the same name)

4	- NCBI species TaxID

5	- species taxonomic lineage

============================================================
===	Supplementary Dataset 2
===	Supplementary_Dataset_2.profiles.tgz
===	Sequence profiles and multiple alignments
============================================================
The gzipped tar archive contains 566 sequence profiles in aligned FASTA format with consensus sequence plus additional multiple alignments in two subdirectories: 21 alignments in Type_V_profiles and 15 alignments in Type_VI_profiles. The latter 36 alignments were used to construct trees in the Supplementary Figures 6 and 7 and are in the simple alignment format (sequence ID, tab or space, sequence). Also, the file profFam.tab (plain ASCII tab-delimited text) contains the description of the sequence profiles.

............................................................
...	profFam.tab
............................................................
1	- profile name
2	- family name
3	- CRISPR-Cas (sub)type

============================================================
===	Supplementary Dataset 3
===	Supplementary_Dataset_3.xlsx
===	Modules derived from bipartite network
============================================================

============================================================
===	Supplementary Dataset 4
===	Supplementary_Dataset_4.trees.tgz
===	Sequence-based trees
============================================================
The gzipped tar archive contains seven files:

Cas1.faa	- multiple alignment of Cas1 sequences in aligned FASTA format
Cas1.tre	- tree of Cas1 sequences in Newick format
Cas1_tree_info.xlsx	- description of Cas1 sequence IDs in Esxcel format
Cas9.tre	- tree of Cas9 sequences in Newick format
Cas9_tree_info.xlsx	- description of Cas9 sequence IDs in Esxcel format
TnpB.tre	- tree of TnpB and Cas12f sequences in Newick format
TnpB_tree_info.xlsx	- description of TnpB and Cas12f sequence IDs in Esxcel format

============================================================
===	Supplementary Dataset 5
===	Supplementary_Dataset_5.xlsx
===	CRISPR-Cas loci corresponding to new subtypes and yet unclassified systems
============================================================

============================================================
===	Supplementary Dataset 6
===	Supplementary_Dataset_6.distr.tgz
===	Distribution of CRISPR-Cas loci across prokaryotic taxa
============================================================
The gzipped tar archive contains two plain ASCII tab-delimited text files:

cas1903.tax.nn.txt	- number of CRISPR-Cas loci of a given subtype in a given taxonomic group

cas1903.tax.wt.txt	- sum of genome assembly weights for CRISPR-Cas loci of a given subtype in a given taxonomic group
