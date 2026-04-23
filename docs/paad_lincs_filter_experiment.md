# PAAD LINCS Cell-Filter Experiment

## Question

The colon pipeline filters LINCS L1000 signatures to colon cell lines before creating drug-level LINCS features. For PAAD, we tested the analogous idea under the strict raw-only setting.

The initially available raw LINCS sources had limited PAAD overlap:

| LINCS source | PAAD cell-line overlap |
|---|---:|
| GSE70138 | `YAPC` only |
| GSE92742 | none in the available metadata |

Therefore the PAAD-specific LINCS test is a `YAPC`-only sensitivity experiment, not a broad PAAD cell-line panel like colon.

We then added CMap/LINCS2020 beta as an external raw LINCS source. It still adds only `YAPC` as a direct PAAD/pancreatic overlap cell line, but it provides many more YAPC perturbation signatures.

## Runs

| Condition | Config | Run ID |
|---|---|---|
| All-cell raw LINCS | `configs/paad_raw_only.json` | `260423_PAAD_V2_raw_only` |
| GSE70138 PAAD-overlap LINCS | `configs/paad_raw_lincs_yapc_only.json` | `260423_PAAD_V2_raw_lincs_yapc_only` |
| CMap2020 PAAD-overlap LINCS | `configs/paad_raw_lincs_cmap2020_yapc.json` | `260423_PAAD_V2_raw_lincs_cmap2020_yapc` |
| No LINCS control | `configs/paad_raw_no_lincs.json` | `260423_PAAD_V2_raw_no_lincs` |

All runs use raw GDSC2, raw DepMap, raw DrugBank/LINCS SMILES matching, and the same 4-fold random/drug/scaffold benchmark.

## LINCS Coverage

| LINCS condition | Cell filter | LINCS-mapped drugs | Raw signatures used |
|---|---|---:|---:|
| all-cell LINCS | all matched `trt_cp` signatures | 120 | 17,978 |
| GSE70138 YAPC-only LINCS | `cell_id == YAPC` | 85 | 1,147 |
| CMap2020 YAPC-only LINCS | `cell_iname == YAPC` | 116 | 2,568 |
| no LINCS | none | 0 | 0 |

## Best Split Results

| LINCS condition | Best random Spearman | Best drug-group Spearman | Best scaffold Spearman |
|---|---:|---:|---:|
| all-cell LINCS | 0.9127 | 0.7123 | 0.7061 |
| GSE70138 YAPC-only LINCS | 0.9110 | 0.6410 | 0.6298 |
| CMap2020 YAPC-only LINCS | 0.9052 | 0.6887 | 0.6530 |
| no LINCS | 0.9107 | 0.4989 | 0.4634 |

Detailed comparison report:

```text
outputs/reports/260423_PAAD_V2_lincs_filter_comparison/
```

S3:

```text
s3://say2-4team/20260409_eseo/260423_PAAD_V2/lincs_filter_comparison/reports/
```

## Interpretation

CMap2020 gives a stronger biology-aligned PAAD-overlap LINCS source than GSE70138 YAPC-only: 116 matched drugs and 2,568 matched signatures after GDSC drug matching. This improves drug-group and scaffold performance relative to the first YAPC-only test.

The all-cell LINCS aggregation remains the strongest strict raw-only setting for PAAD. The YAPC-only filters are biologically stricter and match the colon-style logic, but PAAD-specific LINCS remains a one-cell-line sensitivity experiment because both GSE70138 and CMap2020 only overlap directly on YAPC.

The no-LINCS control drops sharply on drug-group and scaffold splits. This suggests LINCS contributes meaningful generalization signal, but the best PAAD setting is currently the raw-derived all-cell drug-level LINCS aggregation, not the YAPC-only cancer-filtered variant.

For reporting:

```text
LINCS features were generated from raw LINCS Level5 data as drug-level perturbation signatures. PAAD-overlap sensitivity tests restricted LINCS to YAPC, the only direct pancreatic overlap cell line found in GSE70138 and CMap/LINCS2020. CMap2020 YAPC-only improved over GSE70138 YAPC-only, but all-cell drug-level LINCS aggregation remained the strongest overall setting, so the final PAAD raw-only benchmark uses all-cell LINCS while documenting YAPC-only controls.
```
