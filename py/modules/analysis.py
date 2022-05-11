from statistics import mean
from typing import Dict, List, Tuple
from py.models.GroundTruth import Genome, GroundTruth
from py.models.GenePrediction import GenePredictionInfo, GenePredictionResults
from py.constants import MIN_DOMAIN_TOLERANCE, MAX_DOMAIN_TOLERANCE
from py.modules.profile_map import parse_profile_family_map


class CasFamilyCount:
    actual_count: int
    true_positives: int
    false_positives: int

    def __init__(self) -> None:
        self.actual_count = 0
        self.true_positives = 0
        self.false_positives = 0

    def false_negatives(self) -> int:
        return self.actual_count - self.true_positives


def precision(TP: float, FP: float) -> float:
    if TP == 0:
        return 0
    return TP / (TP + FP)


def recall(TP: float, FN: float) -> float:
    if TP == 0:
        return 0
    return TP / (TP + FN)


def accuracy(TP: float, FP: float, FN: float) -> float:
    if TP == 0:
        return 0
    return TP / (TP + FP + FN)


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
) -> Tuple[float, float, List[GenePredictionInfo], List[GenePredictionInfo]]:
    """
    For a given genome and a set of predictions about that genome,
    find the precision and recall of the prediction results
    compared to the ground truth about that genome.

    returns (precision, recall, TPs)
    """

    # Gene IS in groundtruth AND IS predicted
    true_positives: List[GenePredictionInfo] = []
    # Gene IS NOT in groundtruth AND IS predicted
    false_positives: List[GenePredictionInfo] = []

    for gene in genome_truth.genes:
        for prediction in predictions:
            if prediction.profile in gene.profiles and not prediction.visited:
                if gene.start_domain == prediction.start_domain and \
                        gene.end_domain == prediction.end_domain:
                    prediction.accuracy = 1.0
                    true_positives.append(prediction)
                else:
                    # Extract start and end domains
                    true_start_domain = gene.start_domain
                    true_end_domain = gene.end_domain

                    predicted_start_domain = prediction.start_domain
                    predicted_end_domain = prediction.end_domain

                    # Calculate domain span (length)
                    true_domain_span = true_end_domain - true_start_domain
                    predicted_domain_span = predicted_end_domain - \
                        predicted_start_domain

                    if predicted_domain_is_within_error_margins(
                            true_start_domain, true_end_domain,
                            predicted_start_domain, predicted_end_domain
                    ):
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

                prediction.visited = True
                break

    false_positives += list(filter(lambda p: not p.visited, predictions))

    num_expected_genes = len(genome_truth.genes)
    num_true_positives = len(true_positives)
    # False negatives are equivalent to the number of genes
    # hmmer DIDN'T predict correctly.
    # E.g. gene IS in groundtruth AND IS NOT predicted
    num_false_negatives = num_expected_genes - num_true_positives
    num_false_positives = len(false_positives)

    p = precision(num_true_positives, num_false_positives)
    r = recall(num_true_positives, num_false_negatives)

    return (p, r, true_positives, false_positives)


def pipeline_statistics(
    groundtruth: GroundTruth,
    prediction_results: GenePredictionResults
) -> Tuple[float, float]:
    """
    Given a set of gene predictions against a ground truth of genomes,
    calculate the average precision and recall of all predictions.
    """

    cas_profile_families = parse_profile_family_map()
    count_of_families: Dict[str, CasFamilyCount] = {}

    precisions = []
    recalls = []
    accuracies = []
    tp_e_vals = []
    fp_e_vals = []

    for genbank_id in prediction_results.results:
        print(f"Analysing results for {genbank_id}...")

        genome = groundtruth.genomes[genbank_id]
        predictions = prediction_results.get_sorted_results(genbank_id)

        (_, _, TPs, FPs) = genome_prediction_statistics(genome, predictions)

        for gene in genome.genes:
            for profile in gene.profiles:
                if profile in cas_profile_families:
                    family = cas_profile_families[profile].family

                    if family not in count_of_families:
                        # Init the key/object pair if it doesn't exist yet
                        count_of_families[family] = CasFamilyCount()

                    count_of_families[family].actual_count += 1

        for tp in TPs:
            family = cas_profile_families[tp.profile].family

            if family not in count_of_families:
                # Init the key/object pair if it doesn't exist yet
                count_of_families[family] = CasFamilyCount()

            count_of_families[family].true_positives += 1

        for fp in FPs:
            family = cas_profile_families[fp.profile].family

            if family not in count_of_families:
                # Init the key/object pair if it doesn't exist yet
                count_of_families[family] = CasFamilyCount()

            count_of_families[family].false_positives += 1

        if len(TPs) > 0:
            # average_tp_e = mean(map(lambda tp: tp.e_val, TPs))
            tp_e_vals += list(map(lambda tp: tp.e_val, TPs))

        if len(FPs) > 0:
            # average_fp_e = mean(map(lambda fp: fp.e_val, FPs))
            fp_e_vals += list(map(lambda fp: fp.e_val, FPs))

    for family in count_of_families:
        (p, r, a) = family_stats(count_of_families[family])
        precisions.append(p)
        recalls.append(r)
        accuracies.append(a)

    print(f"Average E-Value of TP: {mean(tp_e_vals)}")
    print(f"Average E-Value of FP: {mean(fp_e_vals)}")

    write_per_family_statistics_to_file(count_of_families)

    average_precision = mean(precisions) if len(precisions) > 0 else 0.0
    average_recall = mean(recalls) if len(recalls) > 0 else 0.0
    average_accuracy = mean(accuracies) if len(accuracies) > 0 else 0.0

    return (average_precision, average_recall, average_accuracy)


def family_stats(family: CasFamilyCount) -> Tuple[float, float, float]:
    TPs = family.true_positives
    FPs = family.false_positives
    FNs = family.false_negatives()
    p = precision(TPs, FPs)
    r = recall(TPs, FNs)
    a = accuracy(TPs, FPs, FNs)

    return (p, r, a)


def write_per_family_statistics_to_file(
    family_counts: Dict[str, CasFamilyCount]
) -> None:
    table = "Family,Precision,Recall,Accuracy\n"

    for family in family_counts:
        (p, r, a) = family_stats(family_counts[family])
        table += f"{family},{p},{r},{a}\n"

    f = open("family_statistics.csv", 'r+')
    f.truncate(0)
    f.write(table)
    f.close()
