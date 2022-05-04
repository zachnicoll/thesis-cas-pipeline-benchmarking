import csv
import json
from os.path import exists
import re

from models.GroundTruth import Gene, GroundTruth
from constants import (
    GROUNDTRUTH_FILENAME,
    GROUNDTRUTH_JSON_FILENAME,
    LOCI_HEADER_GENBANK_ID_INDEX,
    LOCI_DESCRIPTION_DOMAIN_INDEX,
    LOCI_DESCRIPTION_CRISPR_INDEX,
    LOCI_DESCRIPTION_PROFILES_INDEX,
    LOCI_DESCRIPTION_SEQUENCE_FAMILIES_INDEX,
    LOCI_DESCRIPTION_SYSTEM_SUBTYPE_INDEX,
)


def parse_gene_row(row: "list[str]") -> Gene:
    [start_domain, end_domain] = map(
        int, row[LOCI_DESCRIPTION_DOMAIN_INDEX].split("..")
    )
    is_crispr_cas = row[LOCI_DESCRIPTION_CRISPR_INDEX] == "+"
    profiles = row[LOCI_DESCRIPTION_PROFILES_INDEX].split(",")
    sequence_families = row[LOCI_DESCRIPTION_SEQUENCE_FAMILIES_INDEX].split(
        ",")
    system_subtype = row[LOCI_DESCRIPTION_SYSTEM_SUBTYPE_INDEX]

    return Gene(
        start_domain,
        end_domain,
        is_crispr_cas,
        profiles,
        sequence_families,
        system_subtype,
    )


def parse_groundtruth() -> GroundTruth:
    print("Parsing groundtruth genome data...")

    if exists(GROUNDTRUTH_JSON_FILENAME):
        print(f"Existing {GROUNDTRUTH_FILENAME} file found, loading contents.")

        # Cached JSON file exists, serialize it as a GroundTruth object
        groundtruth_json = json.load(open(GROUNDTRUTH_JSON_FILENAME))
        return GroundTruth.from_json(groundtruth_json)

    groundtruth_text = open(GROUNDTRUTH_FILENAME).read()

    # Split groundtruth text into alternating rows of loci
    # header and description (using Regex groups), where headers
    # start with '==='.
    split_by_genome = re.split("(===.+)\n", groundtruth_text)

    groundtruth = GroundTruth()
    current_genbank_id = None

    for entry in split_by_genome:
        if entry.startswith("==="):
            # Split string by tab to create a table row (list)
            header: list[str] = entry.split("\t")

            # Extract Genbank ID and intialise genome in GroundTruth object
            current_genbank_id = header[LOCI_HEADER_GENBANK_ID_INDEX]
            groundtruth.init_genome(current_genbank_id)
        elif current_genbank_id is not None:
            # Entry is multiple lines containing each gene & its description
            string_rows = entry.splitlines()

            # Split by tab to create table of data
            table = csv.reader(string_rows, delimiter="\t")

            row: list[str]
            for row in table:
                # Parse row and add to GroundTruth object
                gene = parse_gene_row(row)
                groundtruth.add_gene(current_genbank_id, gene)

    # Convert to JSON
    groundtruth_json = json.dumps(groundtruth, default=vars)

    # Dump to file for local caching
    json_file = open(GROUNDTRUTH_JSON_FILENAME, "a")
    json_file.write(groundtruth_json)
    json_file.close()

    return groundtruth
