# PAAD Strict Raw-Only Run

## Correction

`configs/paad_original_core.json` was not strict enough for the requirement "do not use processed inputs." It rebuilt GDSC/DepMap core tables from source files, but still used curated bridge parquet files for SMILES and LINCS.

The strict run is:

```text
configs/paad_raw_only.json
```

This config does not use these processed inputs:

- `drug_features_catalog_brca_reference_20260420.parquet`
- `lincs_drug_signature_pancancer_20260421.parquet`
- `processed_lincs_20260406/*.parquet`
- any prebuilt `labels.parquet`, `drug_features.parquet`, or `drug_target_mapping.parquet`

Derived parquet files are still written under the run directory, but they are generated inside this run from raw sources.

## Raw Inputs

| Source | Raw input |
|---|---|
| GDSC2 response | `s3://say2-4team/oringinal_raw/gdsc/GDSC2-dataset.local.csv` |
| GDSC compounds | `s3://say2-4team/oringinal_raw/gdsc/screened_compounds_rel_8.5.csv` |
| GDSC cell details | `s3://say2-4team/oringinal_raw/gdsc/Cell_Lines_Details.xlsx` |
| DepMap model metadata | `s3://say2-4team/oringinal_raw/depmap/Model.csv` |
| DepMap CRISPR | `s3://say2-4team/oringinal_raw/depmap/CRISPRGeneEffect.csv` |
| DrugBank | `s3://say2-4team/oringinal_raw/drugbank/drugbank_all_full_database.xml.zip` |
| LINCS Level5 | `s3://say2-4team/PAAD_raw/Linc1000/gse70138_raw/GSE70138_Broad_LINCS_Level5_COMPZ_n118050x12328_2017-03-06.gctx.gz` |
| LINCS metadata | raw GSE70138 `cell_info`, `gene_info`, `inst_info`, `pert_info`, `sig_info` txt.gz files |

The bucket prefix is misspelled as `oringinal_raw`; this is the actual location found in S3.

## Raw Processing

The raw-only pipeline builds every modeling table in-run:

1. Filters GDSC2 response rows to TCGA `PAAD`.
2. Maps GDSC Sanger model IDs to DepMap `ModelID` using raw `Model.csv`.
3. Selects PAAD cell-line CRISPR features from raw `CRISPRGeneEffect.csv`.
4. Builds the SMILES catalog directly from:
   - LINCS raw `pert_info` `canonical_smiles`
   - DrugBank XML `calculated-properties` / `experimental-properties` SMILES
   - GDSC raw drug names and synonyms for matching
5. Builds LINCS drug signatures directly from raw GSE70138:
   - maps GDSC drugs to LINCS perturbagens by canonical SMILES/name
   - selects `trt_cp` signatures
   - reads Level5 GCTX expression matrix
   - aggregates landmark gene signatures to drug level
   - keeps the top 768 variable LINCS gene features
6. Creates the same three input policy variants:
   - `legacy_rich_valid_smiles_only`
   - `legacy_rich_all_drugs_zero_smiles`
   - `colonstyle_compact_baseline`

## QC

| Metric | Value |
|---|---:|
| Labels | 7,513 |
| Cell lines | 29 |
| Drugs | 295 |
| Valid SMILES drugs | 160 |
| LINCS-mapped drugs | 120 |
| LINCS raw signatures used | 17,978 |
| CRISPR features | 3,500 |

The valid SMILES count is lower than the staged/original-core run because no curated SMILES bridge is allowed. This is expected under strict raw-only provenance.

## Benchmark Result

Report directory:

```text
outputs/reports/260423_PAAD_V2_raw_only/
```

S3 upload:

```text
s3://say2-4team/20260409_eseo/260423_PAAD_V2/raw_only/
```

Best split winners:

| Split | Best variant | Best model | Spearman | Pearson | RMSE | AUROC |
|---|---|---|---:|---:|---:|---:|
| random_4fold | `legacy_rich_valid_smiles_only` | `lightgbm` | 0.9127 | 0.9356 | 1.0665 | 0.9623 |
| drug_group_4fold | `legacy_rich_valid_smiles_only` | `weighted_top3_ensemble` | 0.7123 | 0.7490 | 2.0122 | 0.8602 |
| scaffold_group_4fold | `legacy_rich_valid_smiles_only` | `weighted_top3_ensemble` | 0.7061 | 0.7406 | 2.0374 | 0.8602 |

Full-coverage comparison:

| Variant | Drug/scaffold mean Spearman | Drug coverage |
|---|---:|---:|
| `legacy_rich_valid_smiles_only` | 0.7092 | 160 / 295 |
| `legacy_rich_all_drugs_zero_smiles` | 0.5844 | 295 / 295 |
| `colonstyle_compact_baseline` | 0.5387 | 295 / 295 |

Interpretation: strict raw-only PAAD still favors the legacy-rich feature family. The highest scores come from valid-SMILES-only, but that variant covers only 160 drugs because no curated SMILES bridge was used. For downstream recommendation over the full 295-drug GDSC set, `legacy_rich_all_drugs_zero_smiles` is the full-coverage raw-only alternative.
