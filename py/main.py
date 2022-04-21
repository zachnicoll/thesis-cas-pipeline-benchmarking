from modules.groundtruth import parse_groundtruth
from pipelines.prodigal_hmmer_pipeline import prodigal_hmmer_pipeline


def main():
    groundtruth = parse_groundtruth()
    prodigal_hmmer_pipeline()


if __name__ == "__main__":
    main()
