class GenePredictionInfo:
    cas_sequence_family: str
    start_domain: int
    end_domain: int
    score: float

    def __init__(self,
                 cas_sequence_family: str,
                 start_domain: int,
                 end_domain: int,
                 score: float) -> None:
        self.cas_sequence_family = cas_sequence_family
        self.start_domain = start_domain
        self.end_domain = end_domain
        self.score = score


class GenePredictionResults:
    results: "dict[str, list[GenePredictionInfo]]"

    def __init__(self) -> None:
        self.results = {}

    def add_result(self, genbank_id: str, gene_info: GenePredictionInfo):
        if genbank_id in self.results:
            self.results[genbank_id].append(gene_info)
        else:
            self.results[genbank_id] = [gene_info]
