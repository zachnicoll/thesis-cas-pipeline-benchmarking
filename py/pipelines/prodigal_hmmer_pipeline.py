from typing import Tuple
from modules.prodigal import run_prodigal
from models.GenePrediction import GenePredictionResults
from modules.hmmsearch import run_hmmsearch

"""
TODO: Account for partial matches of Cas genes
i.e. where gene domain is slightly off but a percentage
of the match is correct. Determine a scoring system for this.
"""


def prodigal_hmmer_pipeline(execute_prodigal=True) -> Tuple[GenePredictionResults, float, float]:
    """
    Runs prodigal to data_acquisition/genomes.fasta into a set of proteins.
    Then, each Cas gene .hmm profile is run against the set of proteins
    using hmmersearch to determine if certain Cas genes exist in the genome.

    Returns (prediction_result, prodigal_run_time, hmmer_run_time)
    """

    prodigal_run_time = run_prodigal() if execute_prodigal else 0
    (prediction_result, hmmer_run_time) = run_hmmsearch()

    return (prediction_result, prodigal_run_time, hmmer_run_time)

    masterlist_name = "data_acquisition/crispr_groundtruth.txt"
    f = open(masterlist_name)
    masterlist = f.read()  # read file as string for str operations

    genome_crispr_sequences = []

    # Number of genes that match CAS_TYPE in this genome
    expected_cas_genes = []

    for genome_name in genome_names:
        # Find the first position in the string where the genome name occurs,
        # this will be the start of that genome's section
        first_occurrence_of_target = masterlist.find(genome_name)

        # The line will start with "=== ", so take 4 off of
        # the index to get to the start of the line
        start_of_line = first_occurrence_of_target - 4

        # Cut all text before the start of the correct genome's section
        section = masterlist[start_of_line: len(masterlist)]

        # Find the index of the start of the next section
        next_section_start = section.find("===", 4)

        # Remove everything after that index so we are left with only this genome
        section = section[0:next_section_start]

        # Convert into array of lines
        section_lines = section.splitlines()

        # Remove first line because it will be a genome header line
        target_section_lines = section_lines[1: len(section_lines)]

        for line in target_section_lines:
            columns = line.split("\t")

            domain = columns[1]
            sequence_family = columns[9]

            if sequence_family != "CRISPR":
                expected_cas_genes.append(sequence_family)

            genome_crispr_sequences.append(
                (domain, sequence_family, genome_name)
            )

    expected_cas_genes = len(set(expected_cas_genes))
    correct_targets = 0  # Number of targets correctly identified as CAS_TYPE genes

    for target in targets:
        # Format the domain like domain_start..domain_end
        domain_str = f"{target[0]}..{target[1]}"

        if (domain_str, target[2], target[3]) in genome_crispr_sequences:
            correct_targets += 1

    # How many targets were correct out of all targets?
    precision = correct_targets / len(targets)

    # Were all the cas genes detected, that should have been detected?
    recall = correct_targets / expected_cas_genes if expected_cas_genes > 0 else 0

    prodigal_time = prodigal_stop - prodigal_start
    hmmer_time = hmmer_stop - hmmer_start

    print(f"Precision: {precision}")
    print(f"Recall: {recall}")
    print(f"prodigal took: {prodigal_time} seconds")
    print(f"hmmsearch took: {hmmer_time} seconds")
    print(f"For a total runtime of: {prodigal_time + hmmer_time} seconds")
