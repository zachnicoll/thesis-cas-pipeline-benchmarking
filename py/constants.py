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
PRODIGAL_INPUT_FILENAME = "data_acquisition/genomes.fasta"
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

# Statistical Analysis
DOMAIN_TOLERANCE = 0.1  # Predicted domains can be within 10% of ground truth
MIN_DOMAIN_TOLERANCE = 1 - DOMAIN_TOLERANCE
MAX_DOMAIN_TOLERANCE = 1 + DOMAIN_TOLERANCE
