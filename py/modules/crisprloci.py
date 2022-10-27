import os
import re
import time
from typing import List, Tuple

from py.constants import CRISPRLOCI_OUTPUT
from py.models.GenePrediction import GenePredictionResults, GenePredictionInfo


def parse_summary_row(row: List[str]) -> Tuple[str, GenePredictionInfo]:
    gene_id = "_".join(row[1].split("_", 2)[:2])
    start_domain: int = int(row[7])
    end_domain: int = int(row[8])
    family: str = row[2]

    prediction = GenePredictionInfo(
        None,
        family,
        start_domain,
        end_domain,
        None,
        None
    )

    return gene_id, prediction


def run_crisprloci() -> Tuple[GenePredictionResults, float]:
    gene_predictions = GenePredictionResults()

    start = time.perf_counter()

    file_path = os.path.realpath(os.path.dirname(__file__))

    os.system(
        f"python CRISPRloci_standalone.py \
        -f {file_path}/data_acquisition/genomes/genomes.fasta \
        -output {file_path}/{CRISPRLOCI_OUTPUT} \
        -st dna \
        -s HMM2019"
    )

    end = time.perf_counter()

    summary_file = open(f"{CRISPRLOCI_OUTPUT}/cas_results.tab", "r")

    for line in summary_file:
        # Convert string to a table row, whitespace delimited
        row = re.split(";", line)

        # Parse row and extract prediction info
        for (id, gene_info) in parse_summary_row(row):
            gene_predictions.add_result(id, gene_info)

    return gene_predictions, start - end
