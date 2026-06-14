#!/usr/bin/env python3

import argparse
from pathlib import Path
import pandas as pd

def load_clinical_data(clinical_file):
    """Load TCGA clinical sample metadata."""

    clinical = pd.read_csv(clinical_file, sep="\t", comment="#")
    required_columns = ["PATIENT_ID", "SAMPLE_TYPE"]

    for col in required_columns:
        if col not in clinical.columns:
            raise ValueError(f"Missing required column '{col}'")

    return clinical

def split_matrix(matrix_file, clinical_df, output_dir):

    print(f"\nLoading: {matrix_file}")

    matrix = pd.read_csv(matrix_file, sep="\t")

    mutation_column = matrix.columns[0]
    sample_columns = list(matrix.columns[1:])

    metadata_samples = set(clinical_df["PATIENT_ID"])

    missing_metadata = sorted(set(sample_columns) - metadata_samples)

    print(f"Samples in matrix: {len(sample_columns)}")
    print(f"Samples missing metadata: {len(missing_metadata)}")

    # if missing_metadata:
    #     print("\nMissing samples:")
    #     for sample in missing_metadata:
    #         print(sample)

    valid_samples = [
        s for s in sample_columns
        if s in metadata_samples]

    sample_type_map = dict(
        zip(
            clinical_df["PATIENT_ID"],
            clinical_df["SAMPLE_TYPE"]))

    primary_samples = [
        s for s in valid_samples
        if sample_type_map[s] == "Primary"]

    metastatic_samples = [
        s for s in valid_samples
        if sample_type_map[s] == "Metastasis"
    ]

    print(f"\nPrimary samples: {len(primary_samples)}")
    print(f"Metastatic samples: {len(metastatic_samples)}")

    primary_matrix = matrix[[mutation_column] + primary_samples]
    metastatic_matrix = matrix[[mutation_column] + metastatic_samples]

    output_dir = Path(output_dir)
    output_dir.mkdir(
        parents=True,
        exist_ok=True
    )

    stem = Path(matrix_file).stem

    primary_file = (output_dir / f"{stem}.primary.all")
    metastatic_file = (output_dir / f"{stem}.metastatic.all")

    primary_matrix.to_csv(primary_file, sep="\t", index=False)

    metastatic_matrix.to_csv(metastatic_file, sep="\t", index=False)

    print("\nSaved:")
    print(primary_file)
    print(metastatic_file)


def main():
    parser = argparse.ArgumentParser(
        description=(
            "Split mutational matrix into Primary and Metastatic cohorts"))

    parser.add_argument(
        "--clinical",
        required=True,
        help="Path to data_clinical_sample.txt")

    parser.add_argument(
        "--input",
        required=True,
        help="Path to *.all matrix")

    parser.add_argument(
        "--output",
        required=True,
        help="Output directory")

    args = parser.parse_args()

    clinical = load_clinical_data(args.clinical)

    split_matrix(args.input, clinical, args.output)


if __name__ == "__main__":
    main()