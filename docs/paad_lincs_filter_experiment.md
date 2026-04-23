# PAAD LINCS Cell-Filter Experiment

## Question

The colon pipeline filters LINCS L1000 signatures to colon cell lines before creating drug-level LINCS features. For PAAD, we tested the analogous idea under the strict raw-only setting.

The available raw LINCS sources have limited PAAD overlap:

| LINCS source | PAAD cell-line overlap |
|---|---:|
| GSE70138 | `YAPC` only |
| GSE92742 | none in the available metadata |

Therefore the PAAD-specific LINCS test is a `YAPC`-only sensitivity experiment, not a broad PAAD cell-line panel like colon.

## Runs

| Condition | Config | Run ID |
|---|---|---|
| All-cell raw LINCS | `configs/paad_raw_only.json` | `260423_PAAD_V2_raw_only` |
| PAAD-overlap LINCS | `configs/paad_raw_lincs_yapc_only.json` | `260423_PAAD_V2_raw_lincs_yapc_only` |
| No LINCS control | `configs/paad_raw_no_lincs.json` | `260423_PAAD_V2_raw_no_lincs` |

All three runs use raw GDSC2, raw DepMap, raw DrugBank/LINCS SMILES matching, and the same 4-fold random/drug/scaffold benchmark.

## LINCS Coverage

| LINCS condition | Cell filter | LINCS-mapped drugs | Raw signatures used |
|---|---|---:|---:|
| all-cell LINCS | all matched `trt_cp` signatures | 120 | 17,978 |
| YAPC-only LINCS | `cell_id == YAPC` | 85 | 1,147 |
| no LINCS | none | 0 | 0 |

## Best Split Results

| LINCS condition | Best random Spearman | Best drug-group Spearman | Best scaffold Spearman |
|---|---:|---:|---:|
| all-cell LINCS | 0.9127 | 0.7123 | 0.7061 |
| YAPC-only LINCS | 0.9110 | 0.6410 | 0.6298 |
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

The all-cell LINCS aggregation remains the strongest strict raw-only setting for PAAD. The YAPC-only filter is biologically stricter and matches the colon-style logic, but it loses coverage because YAPC is the only overlapping PAAD cell line in GSE70138.

The no-LINCS control drops sharply on drug-group and scaffold splits. This suggests LINCS contributes meaningful generalization signal, but the best PAAD setting is currently the raw-derived all-cell drug-level LINCS aggregation, not the YAPC-only cancer-filtered variant.

For reporting:

```text
LINCS features were generated from raw LINCS Level5 data as drug-level perturbation signatures. A PAAD-overlap sensitivity test restricted LINCS to YAPC, the only overlapping PAAD LINCS cell line available in GSE70138. The YAPC-only setting improved over no-LINCS but underperformed all-cell drug-level LINCS aggregation, so the final PAAD raw-only benchmark uses all-cell LINCS while documenting the YAPC-only control.
```
