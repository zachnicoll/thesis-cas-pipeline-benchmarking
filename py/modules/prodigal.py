import os
import time
from py.constants import GENOME_INPUT, PRODIGAL_OUTPUT_FILENAME


def run_prodigal() -> float:
    """
    Runs the prodigal program to find all proteins
    in the data_acquisition/genomes/genomes.fasta file. Result
    is stored in proteins.faa, to be used with hmmer.

    Returns run time of prodigal.
    """

    prodigal_start = time.perf_counter()

    os.system(
        f"prodigal \
        -i {GENOME_INPUT} \
        -a {PRODIGAL_OUTPUT_FILENAME} \
        > /dev/null"  # Silence standard output
    )

    prodigal_stop = time.perf_counter()

    return prodigal_stop - prodigal_start
