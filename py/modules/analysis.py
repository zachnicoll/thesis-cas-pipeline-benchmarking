from statistics import mean
from typing import Dict, List, Tuple
from py.models.CasFamilyCount import CasFamilyCount
from py.models.CasProfileFamily import CasProfileFamilyMap
from py.models.GroundTruth import Genome, GroundTruth, Gene
from py.models.GenePrediction import GenePredictionInfo, GenePredictionResults
from py.constants import MIN_DOMAIN_TOLERANCE, MAX_DOMAIN_TOLERANCE
from py.modules.profile_map import parse_profile_family_map


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

    return (
            min_acceptable_start <= pred_start <= max_acceptable_start and
            min_acceptable_end <= pred_end <= max_acceptable_end
    )


def gene_tp_fp(
        gene: Gene,
        pred: GenePredictionInfo
) -> Tuple[List[GenePredictionInfo], List[GenePredictionInfo]]:
    # Gene IS in ground truth AND IS predicted
    true_positives: List[GenePredictionInfo] = []
    # Gene IS NOT in ground truth AND IS predicted
    false_positives: List[GenePredictionInfo] = []

    # Ensure that prediction is in the same domain as the gene
    if gene.start_domain == pred.start_domain and \
            gene.end_domain == pred.end_domain:
        pred.accuracy = 1.0
        true_positives.append(pred)
    else:
        # Extract start and end domains
        true_start_domain = gene.start_domain
        true_end_domain = gene.end_domain

        predicted_start_domain = pred.start_domain
        predicted_end_domain = pred.end_domain

        # Calculate domain span (length)
        true_domain_span = true_end_domain - true_start_domain
        predicted_domain_span = predicted_end_domain - predicted_start_domain

        if predicted_domain_is_within_error_margins(
                true_start_domain, true_end_domain,
                predicted_start_domain, predicted_end_domain
        ):
            # Ensure that accuracy is always <= 1.0 by selecting
            # the larger denominator
            acc = (
                predicted_domain_span / true_domain_span
                if true_domain_span >= predicted_domain_span
                else true_domain_span / predicted_domain_span
            )
            pred.accuracy = acc

            # Prediction was correct, within error margins
            true_positives.append(pred)
        else:
            # Prediction was not correct or within error margins
            false_positives.append(pred)

    pred.visited = True

    return true_positives, false_positives


def genome_prediction_statistics(
        genome_truth: Genome,
        predictions: List[GenePredictionInfo],
        profile_family_map: CasProfileFamilyMap,
        use_score: bool,
) -> Tuple[List[GenePredictionInfo], List[GenePredictionInfo]]:
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
        if use_score:
            highest_confidence: GenePredictionInfo = None

            # Find the prediction relevant to this gene with the highest confidence score
            for prediction in predictions:
                # Prediction is only correct if the same sequence family is detected
                if (
                        prediction.profile is not None and (
                        profile_family_map[prediction.profile].family in gene.sequence_families) or
                        prediction.family in gene.sequence_families
                ) and not prediction.visited:
                    if (highest_confidence is None or
                            prediction.score > highest_confidence.score):
                        highest_confidence = prediction

            if highest_confidence is not None:
                tps, fps = gene_tp_fp(gene, highest_confidence)

                true_positives += tps
                false_positives += fps

                codon_length = 3

                # Remove all other predictions that fall within this prediction's domain.
                # This reduces false positives, as a given domain can only be classified
                # as a single gene.
                # predictions = list(filter(
                #     lambda p:
                #     not (
                #         # Domain is inside classified domain
                #             (
                #                     p.start_domain >= highest_confidence.start_domain and
                #                     p.end_domain <= highest_confidence.end_domain
                #             ) or
                #             # Domain overlaps the start of classified domain
                #             (
                #                     p.end_domain > highest_confidence.start_domain > p.start_domain
                #             ) or
                #             # Domain overlaps the end of classified domain
                #             (
                #                     p.end_domain > highest_confidence.end_domain > p.start_domain
                #             )
                #     ),
                #     predictions))
        else:
            # Find the prediction relevant to this gene with the highest confidence score
            for prediction in predictions:
                # Prediction is only correct if the same sequence family is detected
                if (
                        profile_family_map[prediction.profile].family in gene.sequence_families or
                        prediction.family in gene.sequence_families
                ) and not prediction.visited:
                    tps, fps = gene_tp_fp(gene, prediction)

                    true_positives += tps
                    false_positives += fps

    false_positives += list(filter(lambda p: not p.visited, predictions))

    return true_positives, false_positives


def pipeline_statistics(
        groundtruth: GroundTruth,
        prediction_results: GenePredictionResults
) -> Tuple[float, float, float]:
    """
    Given a set of gene predictions against a ground truth of genomes,
    calculate the average precision and recall of all predictions.
    """

    cas_profile_families = parse_profile_family_map()
    count_of_families: Dict[str, CasFamilyCount] = {}

    precisions = []
    recalls = []
    accuracies = []

    for genbank_id in prediction_results.results:
        print(f"Analysing results for {genbank_id}...")

        genome = groundtruth.genomes[genbank_id]
        predictions = prediction_results.get_sorted_results(genbank_id)
        use_score = predictions[0].score is not None

        (TPs, FPs) = genome_prediction_statistics(genome, predictions, cas_profile_families, use_score)

        for gene in genome.genes:
            for profile in gene.profiles:
                if profile in cas_profile_families:
                    family = cas_profile_families[profile].family

                    if family not in count_of_families:
                        # Init the key/object pair if it doesn't exist yet
                        count_of_families[family] = CasFamilyCount()

                    count_of_families[family].actual_count += 1

        for tp in TPs:
            family = cas_profile_families[tp.profile].family if tp.family is None else tp.family

            if family not in count_of_families:
                # Init the key/object pair if it doesn't exist yet
                count_of_families[family] = CasFamilyCount()

            count_of_families[family].true_positives += 1

        for fp in FPs:
            family = cas_profile_families[fp.profile].family if fp.family is None else fp.family

            if family not in count_of_families:
                # Init the key/object pair if it doesn't exist yet
                count_of_families[family] = CasFamilyCount()

            count_of_families[family].false_positives += 1

    for family in count_of_families:
        (p, r, a) = family_stats(count_of_families[family])
        precisions.append(p)
        recalls.append(r)
        accuracies.append(a)

    write_per_family_statistics_to_file(count_of_families)

    average_precision = mean(precisions) if len(precisions) > 0 else 0.0
    average_recall = mean(recalls) if len(recalls) > 0 else 0.0
    average_accuracy = mean(accuracies) if len(accuracies) > 0 else 0.0

    return average_precision, average_recall, average_accuracy


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

    f = open("family_statistics.csv", 'w+')
    f.truncate(0)
    f.write(table)
    f.close()
