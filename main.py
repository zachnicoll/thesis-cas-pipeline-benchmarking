import datetime

from py.modules.crisprloci import run_crisprloci
from py.modules.prospector_cas import run_prospector_cas_only
from py.pipelines.prodigal_hmmer_pipeline import prodigal_hmmer_pipeline
from py.modules.groundtruth import parse_groundtruth
from py.modules.analysis import pipeline_statistics

# Min gene length = 80bp
# Max delta between CRISPR and nearest cas gene = 67343bp

# thresholds = [5, 10, 15, 20, 25, 50]
# chunks = [100, 150, 200, 250, 300, 1000]

thresholds = [10]
chunks = [250]

should_run_crisprloci = True
should_run_prospector = False
old_prospector = False
should_run_hmmer = False


def main():
    ground_truth = parse_groundtruth()

    if should_run_crisprloci:
        (crisprloci_predictions, crisprloci_run_time) = run_crisprloci()

        (
            loci_precision,
            loci_recall,
            loci_accuracy
        ) = pipeline_statistics(ground_truth, crisprloci_predictions, f"CRISPRloci")

        print(f"""
                -- CRISPRLloci Pipeline Statistics --
                    Genomes Predicted: {len(crisprloci_predictions.results)}
                    Precision: {loci_precision * 100}%
                    Recall: {loci_recall * 100}%
                    Accuracy: {loci_accuracy * 100}%
                
                    Total Run Time: {crisprloci_run_time}s
                """)

    if should_run_prospector:
        if old_prospector:
            (prospector_predictions, prospector_run_time) = run_prospector_cas_only(0, 0)

            (
                prosp_precision,
                prosp_recall,
                prosp_accuracy
            ) = pipeline_statistics(ground_truth, prospector_predictions, f"old_prosp")

            print(f"""
                    -- [OLD] Prospector Pipeline Statistics --
                        Genomes Predicted: {len(prospector_predictions.results)}
                        Precision: {prosp_precision * 100}%
                        Recall: {prosp_recall * 100}%
                        Accuracy: {prosp_accuracy * 100}%
                    
                        Total Run Time: {prospector_run_time}s
                    """)
        else:
            for c in chunks:
                for t in thresholds:
                    (prospector_predictions, prospector_run_time) = run_prospector_cas_only(t, c)

                    (
                        prosp_precision,
                        prosp_recall,
                        prosp_accuracy
                    ) = pipeline_statistics(ground_truth, prospector_predictions, f"new_prosp_{c}_{t}")

                    print(f"""
                -- [NEW] Prospector Pipeline Statistics --
                    Genomes Predicted: {len(prospector_predictions.results)}
                    Precision: {prosp_precision * 100}%
                    Recall: {prosp_recall * 100}%
                    Accuracy: {prosp_accuracy * 100}%
                
                    Total Run Time: {prospector_run_time}s
                """)

    if should_run_hmmer:
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
        ) = pipeline_statistics(ground_truth, hmmer_predictions)

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
