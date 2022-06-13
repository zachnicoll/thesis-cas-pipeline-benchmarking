
from typing import Dict, List


class CasProfileFamily:
    family: str
    types: List[str]

    def __init__(self, family: str, types: List[str]) -> None:
        self.family = family
        self.types = types


CasProfileFamilyMap = Dict[str, CasProfileFamily]
