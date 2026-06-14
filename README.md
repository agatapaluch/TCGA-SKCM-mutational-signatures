# TCGA-SKCM Mutational Signatures Analysis

This project investigates differences between primary and metastatic melanoma samples from the TCGA-SKCM cohort using mutational signature analysis, machine learning classification, and survival analysis.

## Project overview

The main goal was to evaluate whether mutational signature profiles can distinguish primary from metastatic melanoma and whether these differences are associated with clinical outcomes.

The analysis included:

- division of TCGA-SKCM samples into primary and metastatic groups,
- generation of SBS, DBS, and ID mutational matrices,
- mutational signature assignment using COSMIC reference signatures,
- signature filtering and reassignment using a custom database,
- visualization of mutational signature profiles,
- logistic regression classification of primary vs metastatic samples,
- Kaplan-Meier survival analysis.

## Repository structure

```text
.
├── assignment/              # Signature assignment outputs
├── databases/               # Custom signature databases
├── metadata/                # Clinical/sample metadata
├── plots/                   # Generated plots
├── classifier/              # Machine learning classifier notebook
├── survival_analysis/       # R scripts and survival analysis outputs
├── README.md
└── .gitignore
```
