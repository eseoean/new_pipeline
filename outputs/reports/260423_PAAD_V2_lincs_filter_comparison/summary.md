# PAAD Raw-Only LINCS Source/Filter Comparison

| LINCS condition | Source | Cell filter | LINCS drugs | Matched signatures | Best random Spearman | Best drug-group Spearman | Best scaffold Spearman |
|---|---|---|---|---|---|---|---|
| `all_cell_lincs` | GSE70138 | all | 120 | 17978 | 0.9127 | 0.7123 | 0.7061 |
| `gse70138_yapc_only` | GSE70138 | explicit (YAPC) | 85 | 1147 | 0.9110 | 0.6410 | 0.6298 |
| `cmap2020_yapc_only` | CMap/LINCS2020 | explicit (YAPC) | 116 | 2568 | 0.9052 | 0.6887 | 0.6530 |
| `no_lincs` | none | none | 0 | 0 | 0.9107 | 0.4989 | 0.4634 |

## Split Winners

| LINCS condition | Source | Split | Variant | Model | Spearman | Pearson | RMSE | AUROC |
|---|---|---|---|---|---|---|---|---|
| `all_cell_lincs` | GSE70138 | `drug_group_4fold` | `legacy_rich_valid_smiles_only` | `weighted_top3_ensemble` | 0.7123 | 0.7490 | 2.0122 | 0.8602 |
| `all_cell_lincs` | GSE70138 | `random_4fold` | `legacy_rich_valid_smiles_only` | `lightgbm` | 0.9127 | 0.9356 | 1.0665 | 0.9623 |
| `all_cell_lincs` | GSE70138 | `scaffold_group_4fold` | `legacy_rich_valid_smiles_only` | `weighted_top3_ensemble` | 0.7061 | 0.7406 | 2.0374 | 0.8602 |
| `gse70138_yapc_only` | GSE70138 | `drug_group_4fold` | `legacy_rich_valid_smiles_only` | `weighted_top3_ensemble` | 0.6410 | 0.6955 | 2.1689 | 0.8254 |
| `gse70138_yapc_only` | GSE70138 | `random_4fold` | `legacy_rich_valid_smiles_only` | `lightgbm` | 0.9110 | 0.9337 | 1.0817 | 0.9617 |
| `gse70138_yapc_only` | GSE70138 | `scaffold_group_4fold` | `legacy_rich_valid_smiles_only` | `randomforest` | 0.6298 | 0.6708 | 2.2593 | 0.8225 |
| `cmap2020_yapc_only` | CMap/LINCS2020 | `drug_group_4fold` | `legacy_rich_valid_smiles_only` | `weighted_top3_ensemble` | 0.6887 | 0.7472 | 1.9415 | 0.8549 |
| `cmap2020_yapc_only` | CMap/LINCS2020 | `random_4fold` | `legacy_rich_valid_smiles_only` | `weighted_top3_ensemble` | 0.9052 | 0.9330 | 1.0538 | 0.9615 |
| `cmap2020_yapc_only` | CMap/LINCS2020 | `scaffold_group_4fold` | `legacy_rich_valid_smiles_only` | `weighted_top3_ensemble` | 0.6530 | 0.6877 | 2.1359 | 0.8299 |
| `no_lincs` | none | `drug_group_4fold` | `legacy_rich_valid_smiles_only` | `weighted_top3_ensemble` | 0.4989 | 0.5226 | 2.5751 | 0.7383 |
| `no_lincs` | none | `random_4fold` | `legacy_rich_valid_smiles_only` | `lightgbm` | 0.9107 | 0.9348 | 1.0750 | 0.9614 |
| `no_lincs` | none | `scaffold_group_4fold` | `legacy_rich_all_drugs_zero_smiles` | `weighted_top3_ensemble` | 0.4634 | 0.4911 | 2.3296 | 0.7161 |

## Full-Coverage Policy Check

| LINCS condition | Source | Split | valid-only Spearman | zero-smiles Spearman | valid minus zero |
|---|---|---|---|---|---|
| `all_cell_lincs` | GSE70138 | `drug_group_4fold` | 0.7123 | 0.6058 | +0.1065 |
| `all_cell_lincs` | GSE70138 | `random_4fold` | 0.9127 | 0.7988 | +0.1139 |
| `all_cell_lincs` | GSE70138 | `scaffold_group_4fold` | 0.7061 | 0.5631 | +0.1430 |
| `gse70138_yapc_only` | GSE70138 | `drug_group_4fold` | 0.6410 | 0.5536 | +0.0874 |
| `gse70138_yapc_only` | GSE70138 | `random_4fold` | 0.9110 | 0.7973 | +0.1137 |
| `gse70138_yapc_only` | GSE70138 | `scaffold_group_4fold` | 0.6298 | 0.5443 | +0.0855 |
| `cmap2020_yapc_only` | CMap/LINCS2020 | `drug_group_4fold` | 0.6887 | 0.6100 | +0.0787 |
| `cmap2020_yapc_only` | CMap/LINCS2020 | `random_4fold` | 0.9052 | 0.8133 | +0.0919 |
| `cmap2020_yapc_only` | CMap/LINCS2020 | `scaffold_group_4fold` | 0.6530 | 0.5735 | +0.0794 |
| `no_lincs` | none | `drug_group_4fold` | 0.4989 | 0.4958 | +0.0031 |
| `no_lincs` | none | `random_4fold` | 0.9107 | 0.7960 | +0.1146 |
| `no_lincs` | none | `scaffold_group_4fold` | 0.4590 | 0.4634 | -0.0044 |

## Interpretation

- CMap/LINCS2020 adds a larger PAAD-overlap YAPC source than GSE70138: 116 matched PAAD drugs and 2,568 matched signatures after GDSC drug matching, versus 85 drugs and 1,147 signatures in GSE70138 YAPC-only.
- CMap2020 YAPC-only improves over GSE70138 YAPC-only on drug-group and scaffold splits, but the broad all-cell GSE70138 drug-level LINCS aggregation remains strongest overall on scaffold and slightly stronger on drug-group.
- The no-LINCS control remains much weaker on drug-group and scaffold splits, so LINCS is carrying useful generalization signal.
- Reporting stance: use all-cell raw LINCS as current performance leader, and keep CMap2020 YAPC-only as the stronger biology-aligned PAAD-specific sensitivity run.
