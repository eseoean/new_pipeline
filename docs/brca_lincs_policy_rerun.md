# BRCA LINCS Policy Rerun

Date: 2026-04-23  
Cancer: BRCA, breast invasive carcinoma (유방암)

## Goal

Check whether BRCA is actually an exception to the current practical direction.

BRCA is the one cancer where a strong single-cell precedent already exists:

- historical BRCA protocol: `MCF7`-centered
- current LINCS2020 raw overlap: multiple direct BRCA lines plus deep `MCF7` coverage

So this rerun tested three scopes instead of two:

1. `all-cell`
2. direct BRCA overlap
3. explicit `MCF7-only`

## Raw LINCS source used

- Source: LINCS2020 / CMap2020 raw Level5
- Local raw cache:
  - `data/raw_cache/BRCA_raw/external_lincs/LINCS2020/level5_beta_trt_cp_n720216x12328.gctx`
  - `data/raw_cache/BRCA_raw/external_lincs/LINCS2020/cellinfo_beta.txt`
  - `data/raw_cache/BRCA_raw/external_lincs/LINCS2020/siginfo_beta.txt`
  - `data/raw_cache/BRCA_raw/external_lincs/LINCS2020/compoundinfo_beta.txt`
  - `data/raw_cache/BRCA_raw/external_lincs/LINCS2020/geneinfo_beta.txt`

## BRCA label scope

From the current raw BRCA source after pipeline harmonization:

- labels: 13,106
- cell lines: 51
- drugs: 295

## Direct BRCA cell overlap

Direct LINCS2020 overlaps from the current BRCA training cell set were:

- `BT20`
- `BT474`
- `HS578T`
- `MCF7`
- `MDAMB231`
- `MDAMB468`
- `T47D`

This makes BRCA structurally different from the other reruns:

- PAAD: `1` direct overlap
- HNSC: `0`
- LIHC: `2`
- LUNG(NSCLC): `10`
- BRCA: `7` direct overlaps, plus a canonical single representative precedent (`MCF7`)

## Runs

### 1. All-cell run

- Config: `configs/brca_raw_lincs_cmap2020_all.json`
- Run ID: `260423_BRCA_V2_raw_lincs_cmap2020_all`
- Structural policy tag: `all_cell_main`

Coverage:

- labels: 13,106
- BRCA cell lines: 51
- drugs: 295
- valid SMILES drugs: 186
- LINCS mapped drugs: 176
- matched signatures: 81,207
- included LINCS cell lines: 221

Best scores:

- `legacy_rich_valid_smiles_only`
  - random: Spearman `0.8707`
  - drug-group: Spearman `0.6893`
  - scaffold: Spearman `0.7036`

### 2. Direct overlap run

- Config: `configs/brca_raw_lincs_cmap2020_overlap.json`
- Run ID: `260423_BRCA_V2_raw_lincs_cmap2020_overlap`
- Structural policy tag: `cancer_specific_main`

Coverage:

- overlap LINCS cells: `BT20`, `BT474`, `HS578T`, `MCF7`, `MDAMB231`, `MDAMB468`, `T47D`
- LINCS mapped drugs: 152
- matched signatures: 11,340
- included LINCS cell lines: 7

Best scores:

- `legacy_rich_valid_smiles_only`
  - random: Spearman `0.8692`
  - drug-group: Spearman `0.6556`
  - scaffold: Spearman `0.6679`

### 3. MCF7-only explicit run

- Config: `configs/brca_raw_lincs_cmap2020_mcf7.json`
- Run ID: `260423_BRCA_V2_raw_lincs_cmap2020_mcf7`
- Structural note: explicit representative-line benchmark

Coverage:

- overlap LINCS cell: `MCF7`
- LINCS mapped drugs: 148
- matched signatures: 6,224
- included LINCS cell lines: 1

Best scores:

- `legacy_rich_valid_smiles_only`
  - random: Spearman `0.8699`
  - drug-group: Spearman `0.6708`
  - scaffold: Spearman `0.6821`

## Comparison

| Split | All-cell | Direct overlap | MCF7-only |
|---|---:|---:|---:|
| random_4fold | 0.8707 | 0.8692 | 0.8699 |
| drug_group_4fold | 0.6893 | 0.6556 | 0.6708 |
| scaffold_group_4fold | 0.7036 | 0.6679 | 0.6821 |

Mean drug/scaffold Spearman:

- all-cell: `0.6965`
- direct overlap: `0.6617`
- MCF7-only: `0.6764`

## Input policy result

All three runs recommended:

- `legacy_rich_valid_smiles_only`

So BRCA matches the current input-direction result already seen in PAAD, HNSC, liver, and lung.

## Decision

Final practical BRCA setting should be:

- LINCS main: `all-cell`
- input main: `legacy_rich_valid_smiles_only`
- direct BRCA overlap: supporting/sensitivity comparison
- `MCF7-only`: biologically meaningful reference run, but not the best practical main run

## Interpretation

BRCA is the cleanest stress test of the hybrid rule so far.

- Structurally, BRCA has a strong case for disease-specific LINCS.
- Biologically, `MCF7` is a real representative precedent and not an arbitrary single-cell choice.
- Empirically, though, the broader all-cell block still produces the best held-out drug/scaffold performance.

So BRCA does not break the current practical pattern. It actually strengthens the main lesson:

- disease-overlap coverage tells us what to benchmark,
- but final main-setting selection still has to come from the benchmark itself.
