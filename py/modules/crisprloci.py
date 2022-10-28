import os
import re
import time
from typing import List, Tuple

from py.constants import CRISPRLOCI_OUTPUT
from py.models.GenePrediction import GenePredictionResults, GenePredictionInfo

# CRISPRloci names some families differently, so have to compensate with replacements
family_map = {
    "csm5": "csm5gr7",
    "csm4": "csm4gr5",
    "csm3": "csm3gr7",
    "csm2": "csm2gr11",
    "cmr1": "cmr1gr7",
    "cmr3": "cmr3gr5",
    "cmr4": "cmr4gr7",
    "cmr5": "cmr5gr11",
    "cmr6": "cmr6gr7",
    "cse2": "cse2gr11",
    "abieii": "AbiEii",
    "csa5": "csa5gr11",
    "casr": "casR",
    "cas3hd": "cas3HD",
    "csc1": "csc1gr5",
    "csc2": "csc2gr7",
    "cora": "corA",
    "csx10": "csx10gr5"
}

# Some families have been removed from or don't exist in profFam.tsv
removed = ["wyl", "unknown", "cas3-cas2", "csy2", "csy3", "cse1"]


def parse_summary_row(row: List[str]) -> Tuple[str, GenePredictionInfo]:
    gene_id = "_".join(row[1].split("_", 2)[:2])
    start_domain: int = int(row[7])
    end_domain: int = int(row[8])
    family: str = row[2]

    if family in removed:
        return None, None
    elif family in family_map:
        family = family_map[family]

    if re.search('_[0-9]+$', gene_id) is not None:
        # Remove end characters in some genome naming edge cases
        gene_id = gene_id[:len(gene_id) - 4]

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

    file_path = os.path.realpath(os.path.dirname(sys.argv[0]))
    os.system(
        f"python /home/zach/CRISPRloci-1.0.0/CRISPRloci_standalone.py \
        -f {file_path}/data_acquisition/genomes/genomes.fasta \
        -output {file_path}/{CRISPRLOCI_OUTPUT} \
        -st dna \
        -s HMM2019"
    )

    end = time.perf_counter()

    # Collate results
    os.system(f"cat {CRISPRLOCI_OUTPUT}/*/cas_results.tab > {CRISPRLOCI_OUTPUT}/cas_results.tab")
    summary_file = open(f"{CRISPRLOCI_OUTPUT}/cas_results.tab", "r")

    for line in summary_file:
        if len(line) > 1:
            # Convert string to a table row, whitespace delimited
            row = re.split(";", line)

            # Parse row and extract prediction info
            id, gene_info = parse_summary_row(row)

            if id is not None:
                gene_predictions.add_result(id, gene_info)

    return gene_predictions, start - end
