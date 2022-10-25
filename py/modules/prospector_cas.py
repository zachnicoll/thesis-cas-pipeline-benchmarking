import os
import re
import time
from typing import List, Tuple

from py.constants import ROOT_DIR, PROFILE_FAMILY_MAP, CAS_ALIGNMENTS, GENOME_INPUT, PROSPECTOR_OUTPUT, \
    PROSPECTOR_RESULTS
from py.models.GenePrediction import GenePredictionResults, GenePredictionInfo


def parse_summary_row(row: List[str]) -> List[GenePredictionInfo]:
    start_domain: int = int(row[0])
    end_domain: int = int(row[1])
    families: List[str] = list(set(row[2].split(",")))
    score: float = float(row[3])

    predictions: List[GenePredictionInfo] = []
    for family in families:
        if family != "CRISPR":
            predictions.append(GenePredictionInfo(
                family,
                start_domain,
                end_domain,
                score
            ))

    return predictions


def run_prospector_cas_only(proximal_search: bool) -> Tuple[GenePredictionResults, float]:
    """
    """

    gene_predictions = GenePredictionResults()

    prospector_start = time.perf_counter()

    prox_search_str = "--proxSearch" if proximal_search else ""

    os.system(
        f"/home/zach/repos/thesis/prospector/cmake-build-debug/prospector \
        --prof {ROOT_DIR}/{CAS_ALIGNMENTS} \
        --domMap {ROOT_DIR}/{PROFILE_FAMILY_MAP} \
        --genome {ROOT_DIR}/{GENOME_INPUT} \
        --out {ROOT_DIR}/{PROSPECTOR_OUTPUT} \
        --skipSer \
        --casThreshold 5 \
        --casChunkLength 100 \
        {prox_search_str}"
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
                for gene_info in parse_summary_row(row):
                    gene_predictions.add_result(current_genbank_id, gene_info)

    return gene_predictions, prospector_end - prospector_start
