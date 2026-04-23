# Hybrid Pipeline Design

## Why Hybrid?

지금까지 HNSC, liver, lung, colon 계열 실험을 비교하면 두 가지 결론이 같이 나옵니다.

1. Colon-style raw-first 구성은 데이터 출처, QC, 중간 산출물 추적성이 좋습니다.
2. Legacy-rich 모델 입력셋은 LINCS gene-level, CRISPR, target-pair, context feature가 살아 있어 성능이 더 잘 나오는 경우가 많았습니다.

따라서 이 레포의 기본 원칙은 다음입니다.

```text
Raw source-first and QC: colon-style
Model input and benchmark: legacy-rich
```

## Standard Steps

```text
Step 0. Cancer config
Step 1. Raw source inventory and selective download
Step 2. Canonical intermediate tables
Step 3. Legacy-rich feature engineering
Step 4. Input variants
Step 5. Random/drug/scaffold benchmark
Step 6. Variant policy comparison
Step 7. Upload and report
Step 8. External validation, ADMET, KG, ClinicalTrials explanation
```

Step 8 is intentionally separated from model selection. It explains and annotates selected drugs; it should not silently change the model ranking unless a separate post-ranking rescoring experiment is declared.

## Source Tiers

The pipeline can run with two source tiers.

| Tier | Purpose |
|---|---|
| staged | Fast reproduction from already curated raw-folder parquet/CSV files |
| original_core | Rebuild core labels, cell mapping, and CRISPR features from original source files, while documenting any remaining curated bridges |

For PAAD, `configs/paad.json` is the staged reproduction config, and `configs/paad_original_core.json` is the original-core validation config. The original-core run is meant to answer a stricter provenance question: can the model input be reproduced from unprocessed GDSC and DepMap source files, not only from prebuilt parquet tables?

## Canonical Intermediate Tables

The pipeline writes these tables under:

```text
data/processed/runs/{run_id}/step2_intermediate/
```

| Table | Purpose |
|---|---|
| `labels.parquet` | GDSC response labels after cancer filtering |
| `cell_metadata.parquet` | cell line to DepMap mapping and cohort metadata |
| `drug_master.parquet` | GDSC drug annotation plus canonical SMILES |
| `drug_target_mapping.parquet` | target gene symbols parsed from GDSC putative targets |
| `sample_crispr.parquet` | selected DepMap CRISPR gene-effect features |
| `drug_lincs.parquet` | selected LINCS drug-level features |
| `raw_source_manifest.json` | raw source provenance |
| `input_qc.json` | row, drug, cell, mapping, and label QC |

## Feature Blocks

The primary `legacy_rich` feature family includes:

- sample context one-hot features
- DepMap CRISPR dependency features
- drug Morgan fingerprint bits
- RDKit descriptor features
- SMILES TF-IDF SVD features
- LINCS drug-level gene-signature features
- GDSC target/pathway context
- pair-level target dependency features

The compact colon-style control intentionally uses fewer biological/contextual blocks:

- sample CRISPR
- Morgan/RDKit/SMILES features
- five LINCS summary features

## SMILES Policy

SMILES handling is not hard-coded as a universal truth.

| Policy | Meaning |
|---|---|
| `valid_smiles_only` | remove invalid/no-SMILES drugs before model input |
| `all_drugs_zero_smiles` | keep all drugs and set SMILES-derived features to zero for invalid rows |

HNSC and liver showed that this decision can differ by cancer type and by feature family. Every new cancer should compare both policies under identical random/drug/scaffold folds.

## Selection Rule

The default final input is `legacy_rich_valid_smiles_only` unless the benchmark shows a clear drug/scaffold generalization benefit for another variant.

Recommended decision order:

1. reject variants with severe QC warnings;
2. compare `drug_group_4fold` and `scaffold_group_4fold`;
3. use `random_4fold` as a sanity check, not the sole selection criterion;
4. keep compact colon-style as a reproducibility/control baseline even when it is not the top performer.
