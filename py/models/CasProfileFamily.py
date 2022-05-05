
from typing import List


class CasProfileFamily:
    family: str
    types: List[str]

    def __init__(self, family: str, type: List[str]) -> None:
        self.family = family
        self.type = type
