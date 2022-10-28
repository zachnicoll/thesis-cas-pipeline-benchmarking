from datetime import datetime
from pathlib import Path

ROOT_DIR = Path(__file__).resolve().parent.parent

# Groundtruth Files
GROUNDTRUTH_FILENAME = "data_acquisition/crispr_groundtruth.txt"
GROUNDTRUTH_JSON_FILENAME = "groundtruth.json"

# Loci Header Indices
LOCI_HEADER_GENBANK_ID_INDEX = 4

# Loci Description Indices
LOCI_DESCRIPTION_DOMAIN_INDEX = 1
LOCI_DESCRIPTION_CRISPR_INDEX = 7
LOCI_DESCRIPTION_PROFILES_INDEX = 8
LOCI_DESCRIPTION_SEQUENCE_FAMILIES_INDEX = 9
LOCI_DESCRIPTION_SYSTEM_SUBTYPE_INDEX = 10

# Prodigal
GENOME_INPUT = "data_acquisition/genomes"
PRODIGAL_OUTPUT_FILENAME = "proteins.faa"

# hmmer
HMMSEARCH_OUTPUT_FILENAME = "summary.txt"
HMM_DB_FILENAME = "profiles/database.hmm"
HMM_SUMMARY_ID_INDEX = 0
HMM_SUMMARY_E_INDEX = 4
HMM_SUMMARY_SCORE_INDEX = 5
HMM_SUMMARY_DOMAIN_START_INDEX = 19
HMM_SUMMARY_DOMAIN_END_INDEX = 21
# TODO: Determine "good" score threshold
HMM_THRESHOLD_SCORE = 1e-30

# prospector
PROFILE_FAMILY_MAP = "profiles/profFam.tsv"
CAS_ALIGNMENTS = "profiles/alignments"
PROSPECTOR_OUTPUT = "prospector_output"
PROSPECTOR_RESULTS = f"{PROSPECTOR_OUTPUT}/results/"

# CRISPRLoci
CRISPRLOCI_OUTPUT = "crisprloci_output"

# Statistical Analysis
DOMAIN_TOLERANCE = 0.1  # Predicted domains can be within 10% of ground truth
MIN_DOMAIN_TOLERANCE = 1 - DOMAIN_TOLERANCE
MAX_DOMAIN_TOLERANCE = 1 + DOMAIN_TOLERANCE

START_TIME = datetime.now()
