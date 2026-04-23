# PAAD Hybrid Pipeline Report

- Cancer: Pancreatic adenocarcinoma (췌장암)
- Run ID: `260423_PAAD_V2_raw_no_lincs`
- Raw S3: `s3://say2-4team/PAAD_raw`
- Upload S3: `s3://say2-4team/20260409_eseo/260423_PAAD_V2/raw_no_lincs`
- Fold policy: 4-fold random/drug/scaffold

## Input QC

- Labels: 7513 rows, 29 cell lines, 295 drugs
- Valid SMILES drugs: 160 / 295
- DepMap mapped cell lines: 29 / 29
- CRISPR features: 3500
- LINCS mapped drugs: 0

## LINCS Policy

- Strategy: `no_lincs`
- Main-run eligible: no
- Disease-overlap LINCS cells: none
- Representative cells declared: none
- Matched drugs: 0
- Matched signatures: 0
- Reason: LINCS mode is disabled for this run.

## legacy_rich_valid_smiles_only

- Input shape: [4110, 5756]
- Feature count: 4262
- Drugs: 160
- Cell lines: 29
- Invalid SMILES pairs kept: 0

### drug_group_4fold
- weighted_top3_ensemble: Spearman 0.4989, RMSE 2.5751
- extratrees: Spearman 0.4850, RMSE 2.5822
- randomforest: Spearman 0.4837, RMSE 2.5966
- flat_mlp: Spearman 0.4430, RMSE 2.6997
- residual_mlp: Spearman 0.4262, RMSE 2.7126

### random_4fold
- lightgbm: Spearman 0.9107, RMSE 1.0750
- weighted_top3_ensemble: Spearman 0.9089, RMSE 1.1007
- xgboost: Spearman 0.9045, RMSE 1.1231
- randomforest: Spearman 0.8994, RMSE 1.1712
- extratrees: Spearman 0.8961, RMSE 1.2056

### scaffold_group_4fold
- extratrees: Spearman 0.4590, RMSE 2.6303
- weighted_top3_ensemble: Spearman 0.4520, RMSE 2.6459
- randomforest: Spearman 0.4233, RMSE 2.6906
- lightgbm: Spearman 0.4048, RMSE 2.8540
- xgboost: Spearman 0.3889, RMSE 2.8140

## legacy_rich_all_drugs_zero_smiles

- Input shape: [7513, 5756]
- Feature count: 4265
- Drugs: 295
- Cell lines: 29
- Invalid SMILES pairs kept: 3403

### drug_group_4fold
- weighted_top3_ensemble: Spearman 0.4958, RMSE 2.2541
- randomforest: Spearman 0.4798, RMSE 2.2763
- extratrees: Spearman 0.4708, RMSE 2.2978
- xgboost: Spearman 0.4625, RMSE 2.3148
- residual_mlp: Spearman 0.4467, RMSE 2.3138

### random_4fold
- lightgbm: Spearman 0.7960, RMSE 1.3413
- weighted_top3_ensemble: Spearman 0.7813, RMSE 1.3948
- xgboost: Spearman 0.7811, RMSE 1.4027
- randomforest: Spearman 0.7213, RMSE 1.5481
- flat_mlp: Spearman 0.7074, RMSE 1.5928

### scaffold_group_4fold
- weighted_top3_ensemble: Spearman 0.4634, RMSE 2.3296
- randomforest: Spearman 0.4290, RMSE 2.3670
- xgboost: Spearman 0.4225, RMSE 2.3735
- flat_mlp: Spearman 0.4030, RMSE 2.4416
- lightgbm: Spearman 0.3987, RMSE 2.4238

## colonstyle_compact_baseline

- Input shape: [7513, 5756]
- Feature count: 4678
- Drugs: 295
- Cell lines: 29
- Invalid SMILES pairs kept: 3403

### drug_group_4fold
- weighted_top3_ensemble: Spearman 0.4755, RMSE 2.2756
- extratrees: Spearman 0.4634, RMSE 2.3204
- randomforest: Spearman 0.4606, RMSE 2.3142
- flat_mlp: Spearman 0.4169, RMSE 2.3262
- xgboost: Spearman 0.4145, RMSE 2.3961

### random_4fold
- lightgbm: Spearman 0.7542, RMSE 1.4697
- xgboost: Spearman 0.7416, RMSE 1.5119
- weighted_top3_ensemble: Spearman 0.7405, RMSE 1.5069
- randomforest: Spearman 0.7013, RMSE 1.6128
- flat_mlp: Spearman 0.6999, RMSE 1.6159

### scaffold_group_4fold
- weighted_top3_ensemble: Spearman 0.4352, RMSE 2.3771
- extratrees: Spearman 0.4202, RMSE 2.4023
- randomforest: Spearman 0.4119, RMSE 2.4144
- xgboost: Spearman 0.4099, RMSE 2.4187
- residual_mlp: Spearman 0.3924, RMSE 2.5105

## Variant Policy Comparison

| Split | Best valid-only | Best all-drugs-zero | Delta valid-zero |
|---|---:|---:|---:|
| drug_group_4fold | 0.4989 | 0.4958 | +0.0031 |
| random_4fold | 0.9107 | 0.7960 | +0.1146 |
| scaffold_group_4fold | 0.4590 | 0.4634 | -0.0044 |

## Recommended Input Policy

- Recommended variant: `legacy_rich_all_drugs_zero_smiles`
- Basis: highest mean Spearman across drug_group_4fold and scaffold_group_4fold

| Variant | Mean drug/scaffold Spearman |
|---|---:|
| legacy_rich_all_drugs_zero_smiles | 0.4796 |
| legacy_rich_valid_smiles_only | 0.4789 |
| colonstyle_compact_baseline | 0.4553 |

## Extended Metrics For Split Winners

| Split | Variant | Model | Spearman | Pearson | Kendall | RMSE | MAE | R2 | AUROC | AUPRC |
|---|---|---|---:|---:|---:|---:|---:|---:|---:|---:|
| drug_group_4fold | legacy_rich_valid_smiles_only | weighted_top3_ensemble | 0.4989 | 0.5226 | 0.3527 | 2.5751 | 1.9776 | 0.2712 | 0.7383 | 0.6397 |
| random_4fold | legacy_rich_valid_smiles_only | lightgbm | 0.9107 | 0.9348 | 0.7502 | 1.0750 | 0.8300 | 0.8730 | 0.9614 | 0.9408 |
| scaffold_group_4fold | legacy_rich_all_drugs_zero_smiles | weighted_top3_ensemble | 0.4634 | 0.4911 | 0.3233 | 2.3296 | 1.7772 | 0.2376 | 0.7161 | 0.5407 |

Binary metrics use `label_binary=1` as sensitive and score `-predicted LN_IC50`.
