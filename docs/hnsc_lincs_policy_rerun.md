# HNSC LINCS Policy Rerun

Date: 2026-04-23  
Cancer: HNSC, head and neck squamous cell carcinoma (두경부 편평세포암)

## Goal

Re-run HNSC in `new_pipeline` under the adaptive LINCS policy that was added after PAAD:

- `cancer_specific_main` if overlapping disease LINCS cell lines are sufficient
- `single_representative_main` only when a single canonical representative line has deep coverage
- otherwise `all_cell_main`

## Raw LINCS source used

- Source: LINCS2020 / CMap2020 raw Level5
- Local raw cache:
  - `data/raw_cache/HNSC_raw/external_lincs/LINCS2020/level5_beta_trt_cp_n720216x12328.gctx`
  - `data/raw_cache/HNSC_raw/external_lincs/LINCS2020/cellinfo_beta.txt`
  - `data/raw_cache/HNSC_raw/external_lincs/LINCS2020/siginfo_beta.txt`
  - `data/raw_cache/HNSC_raw/external_lincs/LINCS2020/compoundinfo_beta.txt`
  - `data/raw_cache/HNSC_raw/external_lincs/LINCS2020/geneinfo_beta.txt`

## Runs

### 1. All-cell main candidate

- Config: `configs/hnsc_raw_lincs_cmap2020_all.json`
- Run ID: `260423_HNSC_V2_raw_lincs_cmap2020_all`
- Policy result: `all_cell_main`

Coverage:

- labels: 9,358
- HNSC cell lines: 39
- drugs: 295
- valid SMILES drugs: 186
- LINCS mapped drugs: 176
- matched signatures: 81,207
- included LINCS cell lines: 221
- direct disease-overlap LINCS cells: 0

Best scores:

- `legacy_rich_valid_smiles_only`
  - random: Spearman `0.9091`
  - drug-group: Spearman `0.7317`
  - scaffold: Spearman `0.7063`

### 2. BICR6 explicit supporting run

- Config: `configs/hnsc_raw_lincs_cmap2020_bicr6.json`
- Run ID: `260423_HNSC_V2_raw_lincs_cmap2020_bicr6`
- Policy result: `cancer_specific_unavailable`

Why unavailable:

- CMap2020 contains `BICR6` as a head-and-neck-like line
- but the current HNSC GDSC training cell set does not directly overlap with `BICR6`
- therefore the pipeline does not treat this as a true disease-overlap main-run candidate

Coverage:

- LINCS mapped drugs: 38
- matched signatures: 114
- included LINCS cell lines: 1 (`BICR6`)

Best scores:

- `legacy_rich_valid_smiles_only`
  - random: Spearman `0.9056`
  - drug-group: Spearman `0.5720`
  - scaffold: Spearman `0.5765`

## Comparison

Comparing best split winners:

| Split | All-cell | BICR6 explicit | Delta |
|---|---:|---:|---:|
| random_4fold | 0.9091 | 0.9056 | +0.0035 |
| drug_group_4fold | 0.7317 | 0.5720 | +0.1597 |
| scaffold_group_4fold | 0.7063 | 0.5765 | +0.1297 |

## Input policy result

For both runs, the best input variant was:

- `legacy_rich_valid_smiles_only`

This is important because it differs from the earlier HNSC result where keeping zero-SMILES drugs was better. Under the current raw CMap2020 LINCS rerun, dropping invalid/no-SMILES drugs is clearly stronger than zero-filling them.

## Decision

HNSC should currently use:

- main LINCS policy: `all_cell_main`
- main input policy: `legacy_rich_valid_smiles_only`
- BICR6 explicit: supporting/sensitivity only

Reason:

- no direct GDSC-HNSC and LINCS cell overlap was available
- BICR6-only coverage is too thin to replace the all-cell signal block
- all-cell materially outperforms BICR6 on drug/scaffold holdout
