from modules.groundtruth import parse_groundtruth
from pipelines.prodigal_hmmer_pipeline import prodigal_hmmer_pipeline
from constants import MAX_DOMAIN_TOLERANCE, MIN_DOMAIN_TOLERANCE


def main():
    groundtruth = parse_groundtruth()
    (hmmer_predictions, prodigal_run_time,
     hmmer_run_time) = prodigal_hmmer_pipeline(False)

    for genbank_id in hmmer_predictions.results:
        print(f"Analysing results for {genbank_id}...")

        genome = groundtruth.genomes[genbank_id]
        predictions = hmmer_predictions.get_sorted_results(genbank_id)

        true_positives = []  # Gene IS in groundtruth AND IS predicted
        false_positives = []  # Gene IS NOT in groundtruth AND IS predicted
        true_negatives = []  # Gene IS NOT in groundtruth AND IS NOT predicted
        false_negatives = []  # Gene IS in groundtruth AND IS NOT predicted

        for prediction in predictions:
            gene_candidates = genome.get_candidates_for_prediction(prediction)

            if len(gene_candidates) == 0:
                # There are no genes with the predicted Cas gene family
                # in this genome, therefore prediction is incorrect
                false_positives.append(prediction)
                continue

            for gene in gene_candidates:
                if gene.start_domain == prediction.start_domain and \
                        gene.end_domain == prediction.end_domain:
                    true_positives.append(prediction)
                else:
                    true_start_domain = gene.start_domain
                    true_end_domain = gene.end_domain
                    predicted_start_domain = prediction.start_domain
                    predicted_end_domain = prediction.end_domain

                    true_domain_span = true_end_domain - true_start_domain
                    predicted_domain_span = predicted_end_domain - predicted_start_domain

                    min_acceptable_start = true_start_domain * MIN_DOMAIN_TOLERANCE
                    min_acceptable_end = true_end_domain * MIN_DOMAIN_TOLERANCE

                    max_acceptable_start = true_start_domain * MAX_DOMAIN_TOLERANCE
                    max_acceptable_end = true_end_domain * MAX_DOMAIN_TOLERANCE

                    if predicted_start_domain >= min_acceptable_start and \
                            predicted_end_domain >= min_acceptable_end and \
                            predicted_start_domain <= max_acceptable_start and \
                            predicted_end_domain <= max_acceptable_end:

                        # Predicted gene is inside the allowable domain,
                        # we'll consider it.

                        # Ensure that accuracy is always <= 1.0 by selecting
                        # the larger denominator
                        accuracy = (predicted_domain_span / true_domain_span
                                    if true_domain_span >= predicted_domain_span
                                    else true_domain_span / predicted_domain_span)
                        prediction.accuracy = accuracy

                        # Prediction was correct, within error margins
                        true_positives.append(prediction)
                    else:
                        # Prediction was not correct or within error margins
                        false_positives.append(prediction)

        print(true_positives, false_positives)


if __name__ == "__main__":
    main()
