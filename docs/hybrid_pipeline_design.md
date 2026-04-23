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
| raw_only | Use only unprocessed source files as inputs, then generate every bridge/intermediate table inside this pipeline run |

For PAAD, `configs/paad.json` is the staged reproduction config, `configs/paad_original_core.json` is the original-core validation config, and `configs/paad_raw_only.json` is the strict raw-only config. Raw-only is the correct mode when the rule is "do not use any previously processed bridge parquet as input."

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
| `raw_lincs_build_summary.json` | raw LINCS coverage, overlap cell IDs, matched signatures |
| `lincs_policy_decision.json` | why this cancer uses all-cell, cancer-specific, or supporting-only LINCS |
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

## LINCS Policy

LINCS scope should not be chosen by a single rigid rule such as "single-cell is always bad" or "disease-specific is always better." The policy should prefer the narrowest biologically defensible LINCS scope that still has enough coverage to act as a stable signal block.

The key precedent is BRCA: the reference protocol used `MCF7` only. This is acceptable because `MCF7` is a canonical BRCA line and the LINCS coverage is very deep. By contrast, PAAD currently overlaps mainly on `YAPC`, which is too narrow to replace the stronger all-cell main run.

Recommended decision rule:

1. `cancer_specific_main`
   - Use when two or more overlapping disease LINCS cell lines are available.
2. `single_representative_main`
   - Use when exactly one overlapping disease LINCS cell line is available,
   - and that cell line is a declared representative disease line,
   - and LINCS coverage is still deep enough to be stable.
3. `all_cell_main`
   - Use when disease-specific LINCS overlap is too sparse or too weak.
   - Keep disease-restricted runs as sensitivity analyses.
4. `cancer_specific_supporting_only`
   - Use when a disease-overlap single-cell run exists but is not strong enough to become the main setting.

The pipeline should record this decision explicitly in run outputs with:

- `strategy`
- `main_run_eligible`
- `available_cancer_overlap_cell_ids`
- `matched_drugs`
- `matched_signatures`
- `selection_reason`

Config convention:

```json
"lincs_policy": {
  "representative_cell_lines": ["MCF7"],
  "single_cell_main_min_mapped_drugs": 80,
  "single_cell_main_min_signatures": 5000
}
```

If `representative_cell_lines` is empty, a one-cell disease-specific LINCS run is treated as supporting evidence unless strong justification is added separately.

## Selection Rule

The default final input is `legacy_rich_valid_smiles_only` unless the benchmark shows a clear drug/scaffold generalization benefit for another variant.

Recommended decision order:

1. reject variants with severe QC warnings;
2. compare `drug_group_4fold` and `scaffold_group_4fold`;
3. use `random_4fold` as a sanity check, not the sole selection criterion;
4. keep compact colon-style as a reproducibility/control baseline even when it is not the top performer.
