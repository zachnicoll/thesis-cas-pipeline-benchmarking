from typing import List
from py.models.GenePrediction import GenePredictionInfo


class Gene:
    def __init__(
        self,
        start_domain: int,
        end_domain: int,
        is_crispr_cas: bool,
        profiles: "list[str]",
        sequence_families: "list[str]",
        system_subtype: str,
    ) -> None:
        self.start_domain: int = start_domain
        self.end_domain: int = end_domain
        self.is_crispr_cas: bool = is_crispr_cas
        self.profiles: "list[str]" = profiles
        self.sequence_families: "list[str]" = sequence_families
        self.system_subtype: str = system_subtype

    @staticmethod
    def from_dict(gene: dict) -> 'Gene':
        return Gene(
            gene['start_domain'],
            gene['end_domain'],
            gene['is_crispr_cas'],
            gene['profiles'],
            gene['sequence_families'],
            gene['system_subtype']
        )


class Genome:
    genes: "list[Gene]"

    def __init__(self) -> None:
        self.genes = []

    def get_candidates_for_prediction(
        self,
        prediction: GenePredictionInfo
    ) -> List[Gene]:
        """
        Get all genes in the genome that are of the same Cas sequence
        family.
        """

        relevant_genes = filter(
            lambda gene:
                (prediction.profile in gene.profiles),
            self.genes
        )

        return list(relevant_genes)

    def add_gene(self, gene: Gene) -> None:
        if "CRISPR" not in gene.profiles:
            self.genes.append(gene)


class GroundTruth:
    genomes: "dict[str, Genome]"

    def __init__(self) -> None:
        self.genomes = {}

    def init_genome(self, genbank_id: str) -> None:
        self.genomes[genbank_id] = Genome()

    def add_gene(self, genbank_id: str, gene: Gene) -> None:
        if genbank_id in self.genomes:
            self.genomes[genbank_id].add_gene(gene)
        else:
            raise Exception(
                f"""
                Genbank Id {genbank_id} does not exist yet.
                Did you call init_genome() first?
                """
            )

    @staticmethod
    def from_json(_json: dict) -> 'GroundTruth':
        groundtruth = GroundTruth()
        genomes: dict = _json['genomes']

        genbank_id: str
        for genbank_id in genomes:
            groundtruth.init_genome(genbank_id)
            genome: dict = genomes[genbank_id]

            gene: dict
            for gene in genome["genes"]:
                groundtruth.add_gene(genbank_id, Gene.from_dict(gene))

        return groundtruth
