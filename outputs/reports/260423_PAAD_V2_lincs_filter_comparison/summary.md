# PAAD Raw-Only LINCS Filter Comparison

| LINCS condition | Cell filter | LINCS drugs | LINCS signatures | Best random Spearman | Best drug-group Spearman | Best scaffold Spearman |
|---|---|---:|---:|---:|---:|---:|
| `all_cell_lincs` | `all` | 120 | 17978 | 0.9127 | 0.7123 | 0.7061 |
| `yapc_only_lincs` | `explicit (YAPC)` | 85 | 1147 | 0.9110 | 0.6410 | 0.6298 |
| `no_lincs` | `none` | 0 | 0 | 0.9107 | 0.4989 | 0.4634 |

## Split Winners

| LINCS condition | Split | Variant | Model | Spearman | Pearson | RMSE | AUROC |
|---|---|---|---|---:|---:|---:|---:|
| `all_cell_lincs` | `drug_group_4fold` | `legacy_rich_valid_smiles_only` | `weighted_top3_ensemble` | 0.7123 | 0.7490 | 2.0122 | 0.8602 |
| `all_cell_lincs` | `random_4fold` | `legacy_rich_valid_smiles_only` | `lightgbm` | 0.9127 | 0.9356 | 1.0665 | 0.9623 |
| `all_cell_lincs` | `scaffold_group_4fold` | `legacy_rich_valid_smiles_only` | `weighted_top3_ensemble` | 0.7061 | 0.7406 | 2.0374 | 0.8602 |
| `yapc_only_lincs` | `drug_group_4fold` | `legacy_rich_valid_smiles_only` | `weighted_top3_ensemble` | 0.6410 | 0.6955 | 2.1689 | 0.8254 |
| `yapc_only_lincs` | `random_4fold` | `legacy_rich_valid_smiles_only` | `lightgbm` | 0.9110 | 0.9337 | 1.0817 | 0.9617 |
| `yapc_only_lincs` | `scaffold_group_4fold` | `legacy_rich_valid_smiles_only` | `randomforest` | 0.6298 | 0.6708 | 2.2593 | 0.8225 |
| `no_lincs` | `drug_group_4fold` | `legacy_rich_valid_smiles_only` | `weighted_top3_ensemble` | 0.4989 | 0.5226 | 2.5751 | 0.7383 |
| `no_lincs` | `random_4fold` | `legacy_rich_valid_smiles_only` | `lightgbm` | 0.9107 | 0.9348 | 1.0750 | 0.9614 |
| `no_lincs` | `scaffold_group_4fold` | `legacy_rich_all_drugs_zero_smiles` | `weighted_top3_ensemble` | 0.4634 | 0.4911 | 2.3296 | 0.7161 |

## Full-Coverage Policy Check

| LINCS condition | Split | valid-only Spearman | zero-smiles Spearman | valid minus zero |
|---|---|---:|---:|---:|
| `all_cell_lincs` | `drug_group_4fold` | 0.7123 | 0.6058 | 0.1065 |
| `all_cell_lincs` | `random_4fold` | 0.9127 | 0.7988 | 0.1139 |
| `all_cell_lincs` | `scaffold_group_4fold` | 0.7061 | 0.5631 | 0.1430 |
| `yapc_only_lincs` | `drug_group_4fold` | 0.6410 | 0.5536 | 0.0874 |
| `yapc_only_lincs` | `random_4fold` | 0.9110 | 0.7973 | 0.1137 |
| `yapc_only_lincs` | `scaffold_group_4fold` | 0.6298 | 0.5443 | 0.0855 |
| `no_lincs` | `drug_group_4fold` | 0.4989 | 0.4958 | 0.0031 |
| `no_lincs` | `random_4fold` | 0.9107 | 0.7960 | 0.1146 |
| `no_lincs` | `scaffold_group_4fold` | 0.4590 | 0.4634 | -0.0044 |

## Interpretation

- `all_cell_lincs` is the strongest strict raw-only LINCS setting on drug-group and scaffold generalization.
- `yapc_only_lincs` is the team-style PAAD-overlap filter. It is biologically stricter, but available GSE70138 overlap is only YAPC, reducing LINCS-mapped drugs from 120 to 85 and signatures from 17,978 to 1,147.
- `no_lincs` drops sharply on generalization, especially drug-group and scaffold splits, so LINCS contributes useful signal in this PAAD raw-only pipeline.
- For presentation, describe current best as raw-derived drug-level LINCS all-cell aggregation, and describe YAPC-only as a sensitivity/control experiment rather than the main production setting.
