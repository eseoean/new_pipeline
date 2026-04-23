# Lung LINCS Policy Rerun

Date: 2026-04-23  
Cancer: LUNG, non-small cell lung cancer (비소세포폐암)

## Goal

Re-run lung in `new_pipeline` under the adaptive LINCS policy, with one explicit scope guard:

- use `LUAD + LUSC` only
- exclude `SCLC` completely from the GDSC label universe

This means the current lung rerun is an NSCLC rerun, not a mixed lung-cancer run.

## Raw LINCS source used

- Source: LINCS2020 / CMap2020 raw Level5
- Local raw cache:
  - `data/raw_cache/LUNG_raw/external_lincs/LINCS2020/level5_beta_trt_cp_n720216x12328.gctx`
  - `data/raw_cache/LUNG_raw/external_lincs/LINCS2020/cellinfo_beta.txt`
  - `data/raw_cache/LUNG_raw/external_lincs/LINCS2020/siginfo_beta.txt`
  - `data/raw_cache/LUNG_raw/external_lincs/LINCS2020/compoundinfo_beta.txt`
  - `data/raw_cache/LUNG_raw/external_lincs/LINCS2020/geneinfo_beta.txt`

## NSCLC label scope

From the current raw GDSC source:

- `LUAD`: 15,653 rows
- `LUSC`: 3,863 rows
- combined NSCLC labels: 19,516 raw rows before downstream pair construction
- `SCLC`: excluded

Final rerun QC after pipeline harmonization:

- labels: 19,795
- cell lines: 76
- drugs: 295

## Direct LUAD/LUSC cell overlap

Direct LINCS2020 overlaps from the current LUAD/LUSC training cell set:

- `A549`
- `HCC15`
- `HCC827`
- `NCIH1437`
- `NCIH1563`
- `NCIH1573`
- `NCIH1781`
- `NCIH1975`
- `NCIH838`
- `SKLU1`

So lung is structurally different from:

- PAAD: `1` direct overlap (`YAPC`)
- HNSC: `0` direct overlap
- LIHC: `2` direct overlaps
- LUNG (NSCLC): `10` direct overlaps

This makes lung a strong test case for the hybrid rule:

- structurally, `cancer_specific_main` is clearly available
- empirically, we still need to benchmark it against `all-cell`

## Runs

### 1. All-cell run

- Config: `configs/lung_raw_lincs_cmap2020_all.json`
- Run ID: `260423_LUNG_V2_raw_lincs_cmap2020_all`
- Structural policy tag: `all_cell_main`

Coverage:

- labels: 19,795
- NSCLC cell lines: 76
- drugs: 295
- valid SMILES drugs: 186
- LINCS mapped drugs: 176
- matched signatures: 81,207
- included LINCS cell lines: 221

Best scores:

- `legacy_rich_valid_smiles_only`
  - random: Spearman `0.8978`
  - drug-group: Spearman `0.7430`
  - scaffold: Spearman `0.7147`

### 2. Direct overlap run

- Config: `configs/lung_raw_lincs_cmap2020_overlap.json`
- Run ID: `260423_LUNG_V2_raw_lincs_cmap2020_overlap`
- Structural policy tag: `cancer_specific_main`

Coverage:

- overlap LINCS cells: `A549`, `HCC15`, `HCC827`, `NCIH1437`, `NCIH1563`, `NCIH1573`, `NCIH1781`, `NCIH1975`, `NCIH838`, `SKLU1`
- LINCS mapped drugs: 152
- matched signatures: 5,649
- included LINCS cell lines: 10

Best scores:

- `legacy_rich_valid_smiles_only`
  - random: Spearman `0.8976`
  - drug-group: Spearman `0.7026`
  - scaffold: Spearman `0.7161`

## Comparison

| Split | All-cell | Direct overlap | Delta |
|---|---:|---:|---:|
| random_4fold | 0.8978 | 0.8976 | +0.0002 |
| drug_group_4fold | 0.7430 | 0.7026 | +0.0404 |
| scaffold_group_4fold | 0.7147 | 0.7161 | -0.0015 |

Mean drug/scaffold Spearman:

- all-cell: `0.7288`
- direct overlap: `0.7094`
- delta: `+0.0195`

## Input policy result

For both runs, the best input variant was:

- `legacy_rich_valid_smiles_only`

This means lung follows the same current input direction as PAAD, HNSC, and liver under the hybrid rerun setup.

## Decision

Final practical lung setting should be:

- LINCS main: `all-cell`
- Input main: `legacy_rich_valid_smiles_only`
- direct NSCLC-overlap LINCS: supporting/sensitivity comparison

## Interpretation

Lung is the most interesting case so far:

- unlike PAAD/HNSC, disease-overlap LINCS is not sparse
- unlike liver, the overlap run is not dramatically worse
- but all-cell still wins on the more important aggregate generalization view because drug holdout improves more than the tiny scaffold gain lost by overlap

So the practical lesson holds:

- overlap coverage tells us what disease-specific run is worth testing
- the final main setting still has to be chosen by benchmark, not by structural availability alone
