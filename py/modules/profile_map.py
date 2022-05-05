from typing import Dict
from models.CasProfileFamily import CasProfileFamily


def parse_profile_family_map() -> Dict[str, CasProfileFamily]:
    profile_map_tsv = open("profiles/profFam.tab").read().splitlines()

    # Maps Cas HMM profile name to class containing "family" and "types"
    profile_map: Dict[str, CasProfileFamily] = {}

    for line in profile_map_tsv:
        row = line.split('\t')

        cas_profile = row[0]
        cas_family = row[1]
        cas_types = row[2].split(',')

        profile_map[cas_profile] = CasProfileFamily(cas_family, cas_types)

    return profile_map
