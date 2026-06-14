import pandas as pd
from pathlib import Path

BASE = Path("../assignment/after_filtering")

MATRICES = {
    "SBS": "MutSigMA.SBS96",
    "DBS": "MutSigMA.DBS78",
    "ID": "MutSigMA.ID83"
}

def load_exposures(signature_type, matrix_name):
    """
    Load primary and metastatic assignment activities and return:
        - features dataframe
        - labels series
    """

    primary = pd.read_csv(
        BASE / signature_type / f"{matrix_name}.primary"
        / "Assignment_Solution"
        / "Activities"
        / "Assignment_Solution_Activities.txt",
        sep="\t",
        index_col=0)

    metastatic = pd.read_csv(
        BASE / signature_type / f"{matrix_name}.metastatic"
        / "Assignment_Solution"
        / "Activities"
        / "Assignment_Solution_Activities.txt",
        sep="\t",
        index_col=0)

    labels = pd.concat([
        pd.Series(0, index=primary.index),
        pd.Series(1, index=metastatic.index)])

    df = pd.concat([primary, metastatic], axis=0, sort=False)

    return df, labels

def build_ml_dataset():
    datasets = {}
    labels = None

    for sig_type, matrix_name in MATRICES.items():

        df, current_labels = load_exposures(
            sig_type,
            matrix_name)

        datasets[sig_type] = df

        if labels is None:
            labels = current_labels

    # union of all patients
    all_samples = sorted(
        set().union(
            *[df.index for df in datasets.values()]))

    # align all matrices
    aligned = []

    for df in datasets.values():
        aligned.append(
            df.reindex(all_samples)
        )

    X = pd.concat(
        aligned,
        axis=1
    )

    # missing signature -> 0
    X = X.fillna(0)

    # align labels
    labels = labels.reindex(all_samples)

    if labels.isna().sum() > 0:
        raise ValueError(
            f"Missing labels for {labels.isna().sum()} samples")

    # remove constant features
    constant_cols = X.columns[X.nunique() <= 1]

    X = X.drop(columns=constant_cols)

    print(
        f"Removed {len(constant_cols)} constant features")

    # label last column
    X["label"] = labels.astype(int)

    return X


def main():

    dataset = build_ml_dataset()

    print("\nDataset shape:")
    print(dataset.shape)

    print("\nClass distribution:")
    print(dataset["label"].value_counts())

    print("\nFirst rows:")
    print(dataset.head())

    output_file = "../predictor/ml_dataset.tsv"
    dataset.to_csv(output_file, sep="\t")

    print(
        f"\nSaved dataset to: {output_file}")


if __name__ == "__main__":
    main()