import os
import re
import time
from typing import Tuple
from constants import (
    HMM_DB_FILENAME,
    PRODIGAL_OUTPUT_FILENAME,
    HMMSEARCH_OUTPUT_FILENAME,
    HMM_SUMMARY_ID_INDEX,
    HMM_SUMMARY_DOMAIN_START_INDEX,
    HMM_SUMMARY_DOMAIN_END_INDEX,
    HMM_SUMMARY_SCORE_INDEX,
    HMM_THRESHOLD_SCORE)
from models.GenePrediction import GenePredictionResults, GenePredictionInfo


def parse_summary_row(row: str) \
        -> Tuple[str, GenePredictionInfo]:
    genbank_id: str = re.split(
        r"(.+\.[0-9])", row[HMM_SUMMARY_ID_INDEX])[1]
    start_domain: int = int(row[HMM_SUMMARY_DOMAIN_START_INDEX])
    end_domain: int = int(row[HMM_SUMMARY_DOMAIN_END_INDEX])

    profile: str = row[2]

    # TODO: Use this score to filter out results
    # if score is below threshold
    hmm_score: float = float(row[HMM_SUMMARY_SCORE_INDEX])

    return (
        genbank_id,
        GenePredictionInfo(
            profile,
            start_domain,
            end_domain,
            hmm_score
        )
    )


def run_hmmsearch() -> Tuple[GenePredictionResults, int]:
    print("Executing hmmsearch...")

    gene_predictions = GenePredictionResults()

    hmmer_start = time.perf_counter()

    os.system(
        f"hmmsearch \
            -o /dev/null \
            --tblout {HMMSEARCH_OUTPUT_FILENAME} \
            -E 1.0e-50 \
            {HMM_DB_FILENAME} \
            {PRODIGAL_OUTPUT_FILENAME} \
        "
    )

    summary_file = open(HMMSEARCH_OUTPUT_FILENAME, "r")

    for line in summary_file:
        if line and not line.startswith("#"):
            # Convert string to a table row, whitespace delimited
            row = re.split(r"\s+", line)

            # Parse row and extract prediction info
            (genbank_id, gene_info) = parse_summary_row(row)

            if gene_info.score >= HMM_THRESHOLD_SCORE:
                gene_predictions.add_result(genbank_id, gene_info)

    hmmer_stop = time.perf_counter()
    hmmer_run_time = hmmer_stop - hmmer_start

    return (gene_predictions, hmmer_run_time)
