from py.modules.prospector_cas import run_prospector_cas_only
from py.pipelines.prodigal_hmmer_pipeline import prodigal_hmmer_pipeline
from py.modules.groundtruth import parse_groundtruth
from py.modules.analysis import pipeline_statistics


# Min gene length = 80bp
# Max delta between CRISPR and nearest cas gene = 67343bp

def main():

    groundtruth = parse_groundtruth()

    (prospector_predictions, prospector_run_time) = run_prospector_cas_only(True)

    (
        prosp_precision,
        prosp_recall,
        prosp_accuracy
    ) = pipeline_statistics(groundtruth, prospector_predictions)

    print(f"""
-- Prodigal & hmmer Pipeline Statistics --
    Genomes Predicted: {len(prospector_predictions.results)}
    Precision: {prosp_precision * 100}%
    Recall: {prosp_recall * 100}%
    Accuracy: {prosp_accuracy * 100}%

    Total Run Time: {prospector_run_time}s
""")

    return

    """ Prodigal & hmmer Pipeline """
    execute_prodigal = True  # Only run prodigal when necessary

    (
        hmmer_predictions,
        prodigal_run_time,
        hmmer_run_time
    ) = prodigal_hmmer_pipeline(execute_prodigal)

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

    Prodigal Run Time: {prodigal_run_time}s
    Hmmer Run Time: {hmmer_run_time}s
    Total Runtime: {prodigal_run_time + hmmer_run_time}s
""")


if __name__ == "__main__":
    main()
