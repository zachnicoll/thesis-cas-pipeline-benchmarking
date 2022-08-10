from typing import List, Optional, Dict


class GenePredictionInfo:
    profile: str
    family: Optional[str]
    start_domain: int
    end_domain: int
    score: Optional[float]
    e_val: Optional[float]
    accuracy: Optional[float]
    visited: Optional[bool]

    def __init__(self,
                 profile: str,
                 start_domain: int,
                 end_domain: int,
                 score: float,
                 e_val: float) -> None:
        self.profile = profile
        self.start_domain = start_domain
        self.end_domain = end_domain
        self.score = score
        self.e_val = e_val
        self.accuracy = None
        self.visited = False

    def __init__(self,
                 profile: str,
                 family: str,
                 start_domain: int,
                 end_domain: int) -> None:
        self.profile = profile
        self.family = family
        self.start_domain = start_domain
        self.end_domain = end_domain
        self.score = 0
        self.e_val = 0
        self.accuracy = None
        self.visited = False


class GenePredictionResults:
    results: Dict[str, List[GenePredictionInfo]]

    def __init__(self) -> None:
        self.results = {}

    def add_result(
            self,
            genbank_id: str,
            gene_info: GenePredictionInfo
    ) -> None:
        if genbank_id in self.results:
            # Check if an identical result already exists
            duplicate_result = GenePredictionResults.__duplicate_exists__(
                self.results[genbank_id], gene_info)

            # Only add result if it's unique (hmmer produces duplicates)
            if not duplicate_result:
                self.results[genbank_id].append(gene_info)
        else:
            self.results[genbank_id] = [gene_info]

    def get_sorted_results(self, genbank_id: str) -> List[GenePredictionInfo]:
        """
        Returns the results for a given genome, sorted descending by score.
        """
        self.results[genbank_id].sort(
            reverse=True, key=lambda info: info.score)
        return self.results[genbank_id]

    @staticmethod
    def __compare_results__(
            result_a: GenePredictionInfo,
            result_b: GenePredictionInfo
    ) -> bool:
        """
        Returns true if result_a has the same cas_sequence_family,
        start_domain, and end_domain. Used for detecting duplicate results,
        which hmmer is capable of producing.
        """
        return (result_a.profile == result_b.profile and
                result_a.start_domain == result_b.start_domain and
                result_a.end_domain == result_b.end_domain)

    @staticmethod
    def __duplicate_exists__(
            existing_results: List[GenePredictionInfo],
            new_result: GenePredictionInfo
    ) -> bool:
        for existing_result in existing_results:
            if GenePredictionResults.__compare_results__(
                    existing_result,
                    new_result
            ):
                return True

        return False
