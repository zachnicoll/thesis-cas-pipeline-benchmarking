from pipelines.prodigal_hmmer_pipeline import prodigal_hmmer_pipeline
from modules.groundtruth import parse_groundtruth
from modules.analysis import pipeline_statistics


def main():
    groundtruth = parse_groundtruth()

    """ Prodigal & hmmer Pipeline """
    (
        hmmer_predictions,
        prodigal_run_time,
        hmmer_run_time
    ) = prodigal_hmmer_pipeline(False)

    (
        hmmer_precision,
        hmmer_recall
    ) = pipeline_statistics(groundtruth, hmmer_predictions)

    print(f"""
Prodigal & hmmer pipeline statistics:
    Precision {hmmer_precision * 100}%
    Recall {hmmer_recall * 100}%
    Total Runtime {prodigal_run_time + hmmer_run_time}s
""")


if __name__ == "__main__":
    main()
