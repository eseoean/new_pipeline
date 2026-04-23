# PAAD Original Source Audit

## Summary

The first PAAD V2 run used the existing staged files already present under:

```text
s3://say2-4team/PAAD_raw/
```

That was useful for reproducing the model quickly, but it was not a pure "original raw only" run. The project expectation is stricter: start from unprocessed source files where they exist, then rebuild canonical intermediate tables and model inputs with QC.

For that reason, this repository now includes `configs/paad_original_core.json` and pipeline support for an `original_core` source tier.

Important correction: `original_core` still allows curated bridge tables for SMILES/LINCS. The strict no-processed-input run is documented separately in `docs/paad_raw_only_run.md` and implemented by `configs/paad_raw_only.json`.

## S3 Source Finding

The expected prefix:

```text
s3://say2-4team/original_raw/
```

was empty at inspection time. The actual raw-source prefix in the bucket is misspelled as:

```text
s3://say2-4team/oringinal_raw/
```

The following original source files were copied or referenced from that prefix into PAAD source accounting:

| Source | Original/raw file | PAAD location or handling |
|---|---|---|
| GDSC2 response | `oringinal_raw/gdsc/GDSC2-dataset.local.csv` | copied to `PAAD_raw/original_sources/gdsc/` |
| GDSC cell metadata | `oringinal_raw/gdsc/Cell_Lines_Details.xlsx` | copied to `PAAD_raw/original_sources/gdsc/` |
| GDSC screened compounds | `oringinal_raw/gdsc/screened_compounds_rel_8.5.csv` | copied to `PAAD_raw/original_sources/gdsc/` |
| DepMap model metadata | `oringinal_raw/depmap/Model.csv` | copied to `PAAD_raw/original_sources/depmap/` |
| DepMap CRISPR gene effect | `oringinal_raw/depmap/CRISPRGeneEffect.csv` | copied to `PAAD_raw/original_sources/depmap/` |
| DepMap CRISPR dependency | `oringinal_raw/depmap/CRISPRGeneDependency.csv` | copied to `PAAD_raw/original_sources/depmap/` |
| LINCS metadata | `oringinal_raw/lincs/GSE92742_Broad_LINCS_*_info.txt.gz` | copied to `PAAD_raw/original_sources/lincs/` |
| LINCS Level5 matrix | large GCTX under `oringinal_raw/lincs/` or PAAD LINCS raw | referenced, not duplicated by default |
| ChEMBL SQLite | `oringinal_raw/chembl/chembl_36_sqlite.tar.gz` | referenced, not parsed in the fast core run |
| DrugBank XML | `oringinal_raw/drugbank/full_database.xml` | referenced, not parsed in the fast core run |

## What Original-Core Rebuilds Directly

`configs/paad_original_core.json` rebuilds the PAAD model input using:

- GDSC2 original CSV for drug-cell response labels.
- TCGA `PAAD` filtering from the original GDSC2 table.
- DepMap `Model.csv` for model identity mapping.
- DepMap `CRISPRGeneEffect.csv` for cell-line gene-effect features.
- The same hybrid feature policies used in PAAD V2:
  - `legacy_rich_valid_smiles_only`
  - `legacy_rich_all_drugs_zero_smiles`
  - `colonstyle_compact_baseline`

This run reproduced the PAAD V2 input scale exactly:

| Metric | Value |
|---|---:|
| Labels | 7,513 |
| Cell lines | 29 |
| Drugs | 295 |
| Valid SMILES drugs | 243 |
| LINCS-mapped drugs | 116 |
| CRISPR features | 3,500 |

## Remaining Curated Bridges

Two feature sources are still used as curated bridge tables in `original_core`:

1. Canonical SMILES uses `PAAD_raw/GDSC/drug_features_catalog_brca_reference_20260420.parquet`.
2. LINCS drug-level features use `PAAD_raw/Linc1000/local_enhancements/lincs_drug_signature_pancancer_20260421.parquet`.

The raw files needed to rebuild these bridges are now accounted for or referenced, but parsing them from scratch is intentionally separated from this fast original-core validation because:

- LINCS GCTX is very large.
- ChEMBL and DrugBank name/ID harmonization is a separate source-integration step.
- The immediate purpose here was to verify that PAAD labels, cell mapping, DepMap features, and model-input policy are reproducible from original core sources.

## Original-Core Benchmark Result

The original-core PAAD run produced:

```text
outputs/reports/260423_PAAD_V2_original_core/
```

The recommended input policy remains:

```text
legacy_rich_valid_smiles_only
```

Best split winners from the report:

| Split | Best variant | Best model | Spearman | Pearson | RMSE | AUROC |
|---|---|---|---:|---:|---:|---:|
| random_4fold | `legacy_rich_valid_smiles_only` | `weighted_top3_ensemble` | 0.8887 | 0.9249 | 1.0673 | 0.9550 |
| drug_group_4fold | `legacy_rich_valid_smiles_only` | `weighted_top3_ensemble` | 0.5768 | 0.6142 | 2.1862 | 0.7974 |
| scaffold_group_4fold | `legacy_rich_valid_smiles_only` | `weighted_top3_ensemble` | 0.5758 | 0.5866 | 2.2477 | 0.8001 |

Interpretation: PAAD behaves like liver in this comparison. Keeping only valid-SMILES drugs is better than keeping all drugs with zero-filled SMILES features, and both legacy-rich variants outperform the compact colon-style control on generalization splits.
