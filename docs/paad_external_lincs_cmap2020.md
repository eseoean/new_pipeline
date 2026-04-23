# PAAD External LINCS Search: CMap/LINCS2020

## 목적

PAAD pipeline의 LINCS feature가 기존 GSE70138에서는 `YAPC` 1개 pancreatic cell line에만 겹쳤기 때문에, 외부 LINCS source에서 PAAD/pancreatic cell line coverage를 추가 확보할 수 있는지 확인했다.

## 확보한 원천 source

Source: CMap/LINCS2020 beta public S3 (`macchiato.clue.io/builds/LINCS2020`)

| File | Role |
| --- | --- |
| `cellinfo_beta.txt` | LINCS2020 cell line metadata |
| `siginfo_beta.txt` | signature-level metadata |
| `compoundinfo_beta.txt` | perturbagen/drug metadata, SMILES, target, MOA |
| `geneinfo_beta.txt` | gene metadata and landmark/inferred feature space |
| `level5/level5_beta_trt_cp_n720216x12328.gctx` | compound treatment Level5 expression matrix |

Local cache:

```text
data/raw_cache/PAAD_raw/external_lincs/LINCS2020/
```

Planned S3 raw location:

```text
s3://say2-4team/PAAD_raw/external_lincs/LINCS2020/
```

## PAAD overlap

PAAD GDSC input has 29 cell lines:

```text
AsPC-1, BxPC-3, CAPAN-1, CAPAN-2, CFPAC-1, DAN-G, HPAC, HPAF-II,
Hs-766T, HuP-T3, HuP-T4, KP-1N, KP-2, KP-4, MIA-PaCa-2, MZ1-PC,
PA-TU-8902, PA-TU-8988T, PANC-02-03, PANC-03-27, PANC-04-03,
PANC-08-13, PANC-10-05, PL4, PSN1, SU8686, SUIT-2, SW1990, YAPC
```

CMap/LINCS2020 `cellinfo_beta.txt` and `siginfo_beta.txt` both indicate that the direct PAAD/pancreatic overlap is:

```text
YAPC
```

So CMap2020 does not add new pancreatic cell line names beyond YAPC, but it adds substantially more YAPC drug perturbation signatures than the GSE70138 source used in the first raw-only PAAD run.

## YAPC signature inventory

From `siginfo_beta.txt`, filtered to `pert_type == trt_cp` and `cell_iname == YAPC`:

| Metric | Count |
| --- | ---: |
| YAPC trt_cp signatures | 21,446 |
| Unique CMap perturbagen IDs | 1,787 |
| Unique CMap drug names | 1,760 |
| QC-pass signatures | 12,628 |
| High-quality signatures | 2,590 |
| Exemplar signatures | 1,696 |

Against the current PAAD GDSC drug master:

| Metric | Count |
| --- | ---: |
| PAAD GDSC drugs | 295 |
| Drugs matched to CMap2020 compounds by SMILES/name/alias | 174 |
| Valid SMILES drugs after using CMap2020 + DrugBank raw matching | 186 |
| Drugs with matched YAPC signatures in the pipeline run | 116 |
| Matched YAPC signatures used after GDSC drug matching | 2,568 |

Generated manifests:

```text
data/raw_cache/PAAD_raw/external_lincs/LINCS2020/manifests/paad_cmap2020_cell_overlap.tsv
data/raw_cache/PAAD_raw/external_lincs/LINCS2020/manifests/paad_yapc_sig_manifest.tsv
data/raw_cache/PAAD_raw/external_lincs/LINCS2020/manifests/paad_yapc_compound_summary.tsv
data/raw_cache/PAAD_raw/external_lincs/LINCS2020/paad_cmap2020_overlap_summary.json
```

## Pipeline integration

The hybrid pipeline now supports both older GSE-style LINCS columns and CMap2020 beta columns:

| Concept | GSE70138/GSE92742 style | CMap2020 style |
| --- | --- | --- |
| Drug name | `pert_iname` | `cmap_name` |
| Cell line | `cell_id` | `cell_iname` |
| Gene ID | `pr_gene_id` | `gene_id` |
| Gene symbol | `pr_gene_symbol` | `gene_symbol` |
| Landmark flag | `pr_is_lm` | `feature_space == landmark` |

CMap2020 PAAD-specific run config:

```text
configs/paad_raw_lincs_cmap2020_yapc.json
```

This config keeps the strict raw-only policy and replaces only the LINCS source with CMap/LINCS2020, filtered to `YAPC`.

## Interpretation

CMap2020 is useful as a stronger PAAD-overlap LINCS source because it greatly increases YAPC perturbation coverage. However, it still does not solve the main biological limitation that PAAD-specific LINCS coverage is only one cell line. The fair comparison set is:

1. GSE70138 all-cell drug-level LINCS: broad but not PAAD-specific.
2. GSE70138 YAPC-only LINCS: PAAD-overlap but small.
3. CMap2020 YAPC-only LINCS: PAAD-overlap and much larger.
4. no-LINCS control: measures how much LINCS contributes at all.
