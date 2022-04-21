import os
import time
from constants import PRODIGAL_INPUT_FILENAME, PRODIGAL_OUTPUT_FILENAME


def run_prodigal() -> int:
    """
    Runs the prodigal program to find all proteins
    in the data_acquisition/genomes.fasta file. Result
    is stored in proteins.faa, to be used with hmmer.

    Returns run time of prodigal.
    """

    prodigal_start = time.perf_counter()

    os.system(
        f"prodigal \
        -i {PRODIGAL_INPUT_FILENAME} \
        -a {PRODIGAL_OUTPUT_FILENAME} > /dev/null"
    )

    prodigal_stop = time.perf_counter()

    return prodigal_stop - prodigal_start
