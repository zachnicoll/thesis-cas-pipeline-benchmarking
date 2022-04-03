import os
import pathlib
import re
import time

# 1. Load all Cas gene profiles from profiles folder
profiles = list(pathlib.Path("./profiles").glob("*.hmm"))

# Start timer here to measure speed of prodigal and hmmsearch steps
prodigal_start = time.perf_counter()

# 2. Create list of all proteins in genome with Prodigal
os.system("prodigal -i genomes.fasta -a proteins.faa > /dev/null")

prodigal_stop = time.perf_counter()

# 3. Use hmmer to search the protetins list for Cas proteins

# Store each matched protein in this array as a tuple:
# (target_name, domain_start, domain_end)
targets = []
genome_names = []

print("Running proteins against Cas .hmm profiles...")

hmmer_start = time.perf_counter()

for profile in profiles:
    CAS_TYPE = profile.name.split("_")[0].lower()

    os.system(
        f"hmmsearch --tblout summary.txt {profile.absolute()} proteins.faa > /dev/null"
    )

    # 4. Compare the results with the Cas gene master list to check
    # correctness of matches
    f = open("summary.txt", "r")

    while True:
        line = f.readline()

        if not line:
            break

        if line[0] == "#":
            continue

        # target name is the first column in the summary table, and
        # formatted like "abc.123"
        genome_name = line.split(".")[0]

        if genome_name:
            genome_names.append(genome_name)

        # the domains appear at the end of the line, under the
        # "description of target" column. this regex ensures they are matched
        # quickly and simply.
        match = re.search("# ([0-9]+) # ([0-9]+) #", line)

        if match:
            domain_start = int(match.group(1))
            domain_end = int(match.group(2))

            targets.append((domain_start, domain_end, CAS_TYPE, genome_name))

hmmer_stop = time.perf_counter()

# 5. Check Cas masterlist for correct matches

masterlist_name = "cas1903.isl.type.txt"
f = open(masterlist_name)
masterlist = f.read()  # read file as string for str operations

genome_crispr_sequences = []
expected_cas_genes = []  # Number of genes that match CAS_TYPE in this genome

for genome_name in genome_names:
    # Find the first position in the string where the genome name occurs,
    # this will be the start of that genome's section
    first_occurrence_of_target = masterlist.find(genome_name)

    # The line will start with "=== ", so take 4 off of
    # the index to get to the start of the line
    start_of_line = first_occurrence_of_target - 4

    # Cut all text before the start of the correct genome's section
    section = masterlist[start_of_line : len(masterlist)]

    # Find the index of the start of the next section
    next_section_start = section.find("===", 4)

    # Remove everything after that index so we are left with only this genome
    section = section[0:next_section_start]

    # Convert into array of lines
    section_lines = section.splitlines()

    # Remove first line because it will be a genome header line
    target_section_lines = section_lines[1 : len(section_lines)]

    for line in target_section_lines:
        columns = line.split("\t")

        domain = columns[1]
        sequence_family = columns[9]

        if sequence_family != "CRISPR":
            expected_cas_genes.append(sequence_family)

        genome_crispr_sequences.append((domain, sequence_family, genome_name))

expected_cas_genes = len(set(expected_cas_genes))
correct_targets = 0  # Number of targets correctly identified as CAS_TYPE genes

for target in targets:
    # Format the domain like domain_start..domain_end
    domain_str = f"{target[0]}..{target[1]}"

    if (domain_str, target[2], target[3]) in genome_crispr_sequences:
        correct_targets += 1

# How many targets were correct out of all targets?
precision = correct_targets / len(targets)

# Were all the cas genes detected, that should have been detected?
recall = correct_targets / expected_cas_genes if expected_cas_genes > 0 else 0

prodigal_time = prodigal_stop - prodigal_start
hmmer_time = hmmer_stop - hmmer_start

print(f"Precision: {precision}")
print(f"Recall: {recall}")
print(f"prodigal took: {prodigal_time} seconds")
print(f"hmmsearch took: {hmmer_time} seconds")
print(f"For a total runtime of: {prodigal_time + hmmer_time} seconds")
