from pathlib import Path
import gzip
import re
import pandas as pd

#Ler os datasets e extrair automaticamente os campos disponíveis e transformar em uma tabela
DATASETS = {
    "GSE53697": "data/raw/GSE53697/GSE53697_family.soft.gz",
    "GSE64810": "data/raw/GSE64810/GSE64810_family.soft.gz",
    "GSE68719": "data/raw/GSE68719/GSE68719_family.soft.gz",
}


def parse_soft(path):
    samples = []
    current = None

    with gzip.open(path, "rt", encoding="utf-8", errors="replace") as f:
        for line in f:
            line = line.rstrip("\n")

            if line.startswith("^SAMPLE"):
                if current is not None:
                    samples.append(current)

                gsm = line.split("=", 1)[1].strip()
                current = {
                    "sample_id": gsm,
                    "title": None,
                    "source": None,
                    "characteristics": {},
                }

            elif current is not None and line.startswith("!Sample_title"):
                current["title"] = line.split("=", 1)[1].strip()

            elif current is not None and line.startswith("!Sample_source_name"):
                current["source"] = line.split("=", 1)[1].strip()

            elif current is not None and line.startswith("!Sample_characteristics"):
                value = line.split("=", 1)[1].strip()

                if ":" in value:
                    key, val = value.split(":", 1)
                    current["characteristics"][key.strip().lower()] = val.strip()

        if current is not None:
            samples.append(current)

    rows = []

    for sample in samples:
        row = {
            "sample_id": sample["sample_id"],
            "title": sample["title"],
            "source": sample["source"],
        }

        row.update(sample["characteristics"])
        rows.append(row)

    return pd.DataFrame(rows)


def standardize_metadata(df, study):
    df["study"] = study

    # Criar diagnosis/group de forma específica para cada estudo
    if study == "GSE53697":
        diagnosis = df.get("disease status", pd.Series(index=df.index))

        df["diagnosis"] = diagnosis

        df["group"] = diagnosis.map(
            lambda x: (
                "control"
                if isinstance(x, str) and x.lower() == "control"
                else "disease"
            )
        )

    elif study == "GSE64810":
        diagnosis = df.get("diagnosis", pd.Series(index=df.index))

        df["diagnosis"] = diagnosis

        df["group"] = diagnosis.map(
            lambda x: (
                "control"
                if isinstance(x, str)
                and x.lower() == "neurologically normal"
                else "disease"
            )
        )

    elif study == "GSE68719":
        diagnosis = df.get("diagnosis", pd.Series(index=df.index))

        df["diagnosis"] = diagnosis

        df["group"] = diagnosis.map(
            lambda x: (
                "control"
                if isinstance(x, str)
                and x.lower() == "neurologically normal"
                else "disease"
            )
        )

    return df


def main():
    output_dir = Path("data/processed")
    output_dir.mkdir(parents=True, exist_ok=True)

    for study, soft_path in DATASETS.items():

        print(f"Processing {study}...")

        df = parse_soft(soft_path)
        df = standardize_metadata(df, study)

        study_dir = output_dir / study
        study_dir.mkdir(parents=True, exist_ok=True)

        output = study_dir / "metadata.tsv"

        df.to_csv(output, sep="\t", index=False)

        print(f"  Samples: {len(df)}")
        print(f"  Output: {output}")


if __name__ == "__main__":
    main()