from typing import Tuple
from py.modules.prodigal import run_prodigal
from py.models.GenePrediction import GenePredictionResults
from py.modules.hmmsearch import run_hmmsearch


def prodigal_hmmer_pipeline(
    execute_prodigal=True
) -> Tuple[GenePredictionResults, float, float]:
    """
    Runs prodigal to data_acquisition/genomes.fasta into a set of proteins.
    Then, each Cas gene .hmm profile is run against the set of proteins
    using hmmersearch to determine if certain Cas genes exist in the genome.

    Returns (prediction_result, prodigal_run_time, hmmer_run_time)
    """

    prodigal_run_time = run_prodigal() if execute_prodigal else 0
    (prediction_result, hmmer_run_time) = run_hmmsearch()

    return prediction_result, prodigal_run_time, hmmer_run_time
