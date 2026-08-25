
#produzir uma anotação simplificada e comum aos três estudos.
#com GeneID, Symbol, EnsemblGeneID, GeneType e Description

from pathlib import Path
import gzip
import pandas as pd


DATASETS = [
    "GSE53697",
    "GSE64810",
    "GSE68719",
]


KEEP_COLUMNS = [
    "GeneID",
    "Symbol",
    "EnsemblGeneID",
    "GeneType",
    "Description",
]


def main():

    for study in DATASETS:

        input_file = (
            Path("data/raw")
            / study
            / "Human.GRCh38.p13.annot.tsv.gz"
        )

        output_dir = Path("data/processed") / study
        output_dir.mkdir(parents=True, exist_ok=True)

        output_file = output_dir / "annotation.tsv.gz"

        print(f"Processing {study}...")

        df = pd.read_csv(
            input_file,
            sep="\t",
            compression="gzip",
            dtype=str
        )

        df = df[KEEP_COLUMNS]

        # Remover GeneIDs duplicados
        df = df.drop_duplicates(subset="GeneID")

        df.to_csv(
            output_file,
            sep="\t",
            index=False,
            compression="gzip"
        )

        print(f"  Genes: {len(df)}")
        print(f"  Output: {output_file}")


if __name__ == "__main__":
    main()