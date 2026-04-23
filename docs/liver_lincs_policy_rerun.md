# Liver LINCS Policy Rerun

Date: 2026-04-23  
Cancer: LIHC, liver hepatocellular carcinoma (간암)

## Goal

Check whether liver follows the same final direction as PAAD and HNSC under the adaptive LINCS policy and the current hybrid input design.

## Raw LINCS source used

- Source: LINCS2020 / CMap2020 raw Level5
- Local raw cache:
  - `data/raw_cache/LIHC_raw/external_lincs/LINCS2020/level5_beta_trt_cp_n720216x12328.gctx`
  - `data/raw_cache/LIHC_raw/external_lincs/LINCS2020/cellinfo_beta.txt`
  - `data/raw_cache/LIHC_raw/external_lincs/LINCS2020/siginfo_beta.txt`
  - `data/raw_cache/LIHC_raw/external_lincs/LINCS2020/compoundinfo_beta.txt`
  - `data/raw_cache/LIHC_raw/external_lincs/LINCS2020/geneinfo_beta.txt`

## Direct LIHC cell overlap

From the current GDSC LIHC training set, direct LINCS2020 overlaps were:

- `HUH7`
- `JHH7`

That means liver is structurally different from PAAD and HNSC:

- PAAD: one direct overlap (`YAPC`)
- HNSC: zero direct overlap
- LIHC: two direct overlaps (`HUH7`, `JHH7`)

## Runs

### 1. All-cell run

- Config: `configs/liver_raw_lincs_cmap2020_all.json`
- Run ID: `260423_LIVER_V2_raw_lincs_cmap2020_all`
- Structural policy tag: `all_cell_main`

Coverage:

- labels: 4,164
- LIHC cell lines: 15
- drugs: 295
- valid SMILES drugs: 186
- LINCS mapped drugs: 176
- matched signatures: 81,207
- included LINCS cell lines: 221

Best scores:

- `legacy_rich_valid_smiles_only`
  - random: Spearman `0.8752`
  - drug-group: Spearman `0.6910`
  - scaffold: Spearman `0.6712`

### 2. Direct overlap run

- Config: `configs/liver_raw_lincs_cmap2020_overlap.json`
- Run ID: `260423_LIVER_V2_raw_lincs_cmap2020_overlap`
- Structural policy tag: `cancer_specific_main`

Coverage:

- overlap LINCS cells: `HUH7`, `JHH7`
- LINCS mapped drugs: 61
- matched signatures: 184
- included LINCS cell lines: 2

Best scores:

- `legacy_rich_valid_smiles_only`
  - random: Spearman `0.8680`
  - drug-group: Spearman `0.5590`
  - scaffold: Spearman `0.5435`

## Comparison

| Split | All-cell | Direct overlap | Delta |
|---|---:|---:|---:|
| random_4fold | 0.8752 | 0.8680 | +0.0072 |
| drug_group_4fold | 0.6910 | 0.5590 | +0.1320 |
| scaffold_group_4fold | 0.6712 | 0.5435 | +0.1278 |

## Input policy result

For both runs, the best input variant was:

- `legacy_rich_valid_smiles_only`

This reproduces the earlier liver direction: unlike the earlier HNSC zero-SMILES exception, liver remains stronger when invalid/no-SMILES drugs are removed before model input.

## Decision

Final liver setting should be:

- LINCS main: `all-cell`
- Input main: `legacy_rich_valid_smiles_only`
- Direct-overlap `HUH7 + JHH7`: sensitivity/supporting only

## Interpretation

Liver is an important counterexample for an overly rigid LINCS rule.

- Structurally, the disease-specific run is available and valid.
- Empirically, the all-cell run generalizes much better to held-out drugs/scaffolds.

So the practical lesson is:

- use overlap coverage to decide what disease-specific run is worth testing,
- but choose the final main setting only after direct benchmark comparison.
