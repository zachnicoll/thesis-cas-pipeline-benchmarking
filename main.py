from py.pipelines.prodigal_hmmer_pipeline import prodigal_hmmer_pipeline
from py.modules.groundtruth import parse_groundtruth
from py.modules.analysis import pipeline_statistics


def main():
    groundtruth = parse_groundtruth()

    """ Prodigal & hmmer Pipeline """
    (
        hmmer_predictions,
        prodigal_run_time,
        hmmer_run_time
    ) = prodigal_hmmer_pipeline(False)  # Don't run prodigal while testing

    (
        hmmer_precision,
        hmmer_recall,
        hmmer_accuracy
    ) = pipeline_statistics(groundtruth, hmmer_predictions)

    print(f"""
-- Prodigal & hmmer Pipeline Statistics --
    Genomes Predicted: {len(hmmer_predictions.results)}
    Precision: {hmmer_precision * 100}%
    Recall: {hmmer_recall * 100}%
    Accuracy: {hmmer_accuracy * 100}%
    Total Runtime: {prodigal_run_time + hmmer_run_time}s
""")


if __name__ == "__main__":
    main()
