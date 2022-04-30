from statistics import mean
from typing import List, Tuple
from models.GroundTruth import Genome, GroundTruth
from models.GenePrediction import GenePredictionInfo, GenePredictionResults
from constants import MIN_DOMAIN_TOLERANCE, MAX_DOMAIN_TOLERANCE


def precision(TP: float, FP: float) -> float:
    return TP / (TP + FP)


def recall(TP: float, FN: float) -> float:
    return TP / (TP + FN)


def predicted_domain_is_within_error_margins(
    true_start: int, true_end: int,
    pred_start: int, pred_end: int
) -> bool:
    """
    Determines if the predicted start and end domains of a gene are within
    the allowable error margins for domains.
    See constants.DOMAIN_TOLERANCE.
    """

    # Calculate the acceptable error in predicted domains
    min_acceptable_start = true_start * MIN_DOMAIN_TOLERANCE
    max_acceptable_start = true_start * MAX_DOMAIN_TOLERANCE
    min_acceptable_end = true_end * MIN_DOMAIN_TOLERANCE
    max_acceptable_end = true_end * MAX_DOMAIN_TOLERANCE

    return (pred_start >= min_acceptable_start and
            pred_start <= max_acceptable_start and
            pred_end >= min_acceptable_end and
            pred_end <= max_acceptable_end)


def genome_prediction_statistics(
    genome_truth: Genome,
    predictions: List[GenePredictionInfo]
) -> Tuple[float, float]:
    """
    For a given genome and a set of predictions about that genome,
    find the precision and recall of the prediction results
    compared to the ground truth about that genome.

    returns (precision, recall)
    """

    true_positives = []  # Gene IS in groundtruth AND IS predicted
    false_positives = []  # Gene IS NOT in groundtruth AND IS predicted

    for prediction in predictions:
        gene_candidates = genome_truth.get_candidates_for_prediction(
            prediction)

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
                # Extract start and end domains
                true_start_domain = gene.start_domain
                true_end_domain = gene.end_domain

                predicted_start_domain = prediction.start_domain
                predicted_end_domain = prediction.end_domain

                # Calculate domain span (length)
                true_domain_span = true_end_domain - true_start_domain
                predicted_domain_span = predicted_end_domain - predicted_start_domain

                if predicted_domain_is_within_error_margins(
                        true_start_domain, true_end_domain,
                        predicted_start_domain, predicted_end_domain):
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

    # False negatives are equivalent to the number of genes
    # hmmer DIDN'T predict correctly.
    # E.g. gene IS in groundtruth AND IS NOT predicted
    num_false_negatives = len(genome_truth.genes) - len(true_positives)
    num_true_positives = len(true_positives)
    num_false_positives = len(false_positives)

    p = precision(num_true_positives, num_false_positives)
    r = recall(num_true_positives, num_false_negatives)

    return (p, r)


def pipeline_statistics(
    groundtruth: GroundTruth,
    prediction_results: GenePredictionResults
) -> Tuple[float, float]:
    """
    Given a set of gene predictions against a ground truth of genomes,
    calculate the average precision and recall of all predictions.
    """

    precisions = []
    recalls = []

    for genbank_id in prediction_results.results:
        print(f"Analysing results for {genbank_id}...")

        genome = groundtruth.genomes[genbank_id]
        predictions = prediction_results.get_sorted_results(genbank_id)

        (p, r) = genome_prediction_statistics(genome, predictions)

        precisions.append(p)
        recalls.append(r)

    average_precision = mean(precisions) if len(precisions) > 0 else 0.0
    average_recall = mean(recalls) if len(recalls) > 0 else 0.0

    return (average_precision, average_recall)
