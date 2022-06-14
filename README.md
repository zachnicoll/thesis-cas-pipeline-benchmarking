# Cas Gene Detection Performance Profiling

### Quick Start

Run `download_genomes.R` in RStudio to download a set of genomes. Tweak the `max_genomes` variable to alter the number of genomes downloaded.

Build the HMM database for hmmer by running:

```sh
cd profiles && _build_hmm_db.sh
```

Ensure `execute_prodigal = True` in `main.py` before running the script for the first time.

Then, run a benchmark with:

```py
python main.py
```

## About

This python script carries out the following pipeline for detecting Cas (CRISPR-associated) genes in a given genome:
    
1. Use prodigal to create a list of all the proteins in genome
    
    ```bash
    prodigal -i genome.fna -a proteins.faa
    ```
    
2. Run `hmmsearch` with the created Cas gene profiles, and the list of proteins
    
    ```bash
    hmmsearch --tblout summary.txt profiles/database.hmm proteins.faa
    ```
    
    - The `--tblout` option exports a short summary to the provided filename

4. For each match in the summary file, search the master list of cas genes (`crispr_groundtruth.txt`) for the same domain e.g. `x .. y`. Verify the match by checking if it's the correct type of cas gene (same as type used to generate hmm profile).

5. Collect statistics on how many matches were correctly identified (precision and recall).

### Genome Downloading

`download_genomes.R` downloads a set of genomes in .fasta format from the the NCBI API. All of these genomes are guaranteed to contain some type of Cas gene, as they are pulled straight out of the annotated Cas gene genome file `crispr_groundtruth.txt`.

Doing this allows us to test the script exclusively against genomes that we know *should* contain Cas genes, and so no time is wasted testing the detection flow against genomes that don't contain any Cas genes.

