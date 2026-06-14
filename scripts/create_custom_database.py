#!/usr/bin/env python3

import argparse
import pandas as pd
from pathlib import Path

def detect_relevant_signatures(activities_file, exposure_threshold=0.10, min_patients=5):
    df = pd.read_csv(activities_file, sep="\t")

    sample_col = "Samples"
    sig_cols = [c for c in df.columns if c != sample_col]

    row_sums = df[sig_cols].sum(axis=1)
    exposure = df[sig_cols].div(row_sums, axis=0)
    n_patients = (exposure > exposure_threshold).sum(axis=0)
    keep = n_patients[n_patients >= min_patients].index.tolist()

    return keep


def create_custom_database(cosmic_file, keep_signatures, output_file):

    cosmic = pd.read_csv(cosmic_file, sep="\t")
    first_column = cosmic.columns[0]

    keep_columns = [first_column] + keep_signatures
    missing = [c for c in keep_signatures if c not in cosmic.columns]

    if missing:
        print("\nWarning: signatures not found in COSMIC database: ")
        print(missing)

    keep_columns = [c for c in keep_columns if c in cosmic.columns]
    cosmic_filtered = cosmic[keep_columns]
    cosmic_filtered.to_csv(output_file, sep="\t", index=False)

    return cosmic_filtered.columns.tolist()


def main():

    parser = argparse.ArgumentParser(
        description=(
            "Create custom COSMIC database based on Assignment exposures"))

    parser.add_argument(
        "--activities",
        required=True,
        help="Assignment_Solution_Activities.txt")

    parser.add_argument(
        "--cosmic",
        required=True,
        help="COSMIC reference file")

    parser.add_argument(
        "--output",
        required=True,
        help="Output custom database")

    parser.add_argument(
        "--threshold",
        type=float,
        default=0.10,
        help="Minimum exposure (default=0.10)")

    parser.add_argument(
        "--patients",
        type=int,
        default=5,
        help="Minimum number of patients (default=5)"
    )

    args = parser.parse_args()

    keep_signatures = detect_relevant_signatures(
        args.activities,
        exposure_threshold=args.threshold,
        min_patients=args.patients)

    print(
        f"\nSelected signatures "
        f"({len(keep_signatures)}):"
    )

    for sig in keep_signatures:
        print(sig)

    cols = create_custom_database(
        args.cosmic,
        keep_signatures,
        args.output)

    print(f"\nCustom database saved:")
    print(args.output)
    print(f"\nNumber of columns: {len(cols)-1}")


if __name__ == "__main__":
    main()