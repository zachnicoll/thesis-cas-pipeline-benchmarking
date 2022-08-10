import os
import re
import time
from typing import List, Tuple

from py.constants import ROOT_DIR, PROFILE_FAMILY_MAP, CAS_ALIGNMENTS, GENOME_INPUT, PROSPECTOR_OUTPUT, \
    PROSPECTOR_RESULTS
from py.models.GenePrediction import GenePredictionResults, GenePredictionInfo


def parse_summary_row(row: List[str]) -> GenePredictionInfo:
    start_domain: int = int(row[0])
    end_domain: int = int(row[1])
    families: List[str] = list(set(row[3].split(",")))
    profiles: List[str] = list(set(row[4].split(",")))

    return GenePredictionInfo(
        profiles[0],
        families[0],
        start_domain,
        end_domain
    )


def run_prospector_cas_only() -> Tuple[GenePredictionResults, float]:
    """
    """

    gene_predictions = GenePredictionResults()

    prospector_start = time.perf_counter()

    os.system(
        f"prospector \
        --prof {ROOT_DIR}/{CAS_ALIGNMENTS} \
        --domMap {ROOT_DIR}/{PROFILE_FAMILY_MAP} \
        --genome {ROOT_DIR}/{GENOME_INPUT} \
        --out {ROOT_DIR}/{PROSPECTOR_OUTPUT} \
        --skipSer \
        --casOnly"
    )

    prospector_end = time.perf_counter()

    input_file = os.listdir(GENOME_INPUT)[0].rsplit('.', maxsplit=1)[0]
    summary_file = open(f"{PROSPECTOR_RESULTS}/{input_file}/out_gene.txt", "r")

    current_genbank_id = ""

    for line in summary_file:
        if line and not line.startswith("//"):
            # Convert string to a table row, whitespace delimited
            row = re.split(r"\s+", line)

            if line.startswith(">"):
                current_genbank_id = line.replace(">", "").strip()
            else:
                # Parse row and extract prediction info
                gene_info = parse_summary_row(row)

                gene_predictions.add_result(current_genbank_id, gene_info)

    return gene_predictions, prospector_end - prospector_start
