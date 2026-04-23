# PAAD Hybrid Pipeline Report

- Cancer: Pancreatic adenocarcinoma (췌장암)
- Run ID: `260423_PAAD_V2`
- Raw S3: `s3://say2-4team/PAAD_raw`
- Upload S3: `s3://say2-4team/20260409_eseo/260423_PAAD_V2`
- Fold policy: 4-fold random/drug/scaffold

## Input QC

- Labels: 7513 rows, 29 cell lines, 295 drugs
- Valid SMILES drugs: 243 / 295
- DepMap mapped cell lines: 29 / 29
- CRISPR features: 3500
- LINCS mapped drugs: 116

## legacy_rich_valid_smiles_only

- Input shape: [6219, 6524]
- Feature count: 5036
- Drugs: 243
- Cell lines: 29
- Invalid SMILES pairs kept: 0

### drug_group_4fold
- weighted_top3_ensemble: Spearman 0.5823, RMSE 2.1716
- randomforest: Spearman 0.5493, RMSE 2.2246
- flat_mlp: Spearman 0.5471, RMSE 2.2218
- xgboost: Spearman 0.5457, RMSE 2.2577
- residual_mlp: Spearman 0.5255, RMSE 2.3148

### random_4fold
- lightgbm: Spearman 0.8893, RMSE 1.0549
- weighted_top3_ensemble: Spearman 0.8867, RMSE 1.0799
- xgboost: Spearman 0.8825, RMSE 1.1008
- randomforest: Spearman 0.8728, RMSE 1.1540
- flat_mlp: Spearman 0.8699, RMSE 1.1538

### scaffold_group_4fold
- weighted_top3_ensemble: Spearman 0.5775, RMSE 2.2475
- xgboost: Spearman 0.5654, RMSE 2.2766
- randomforest: Spearman 0.5550, RMSE 2.2895
- lightgbm: Spearman 0.5537, RMSE 2.3136
- flat_mlp: Spearman 0.5321, RMSE 2.3067

## legacy_rich_all_drugs_zero_smiles

- Input shape: [7513, 6524]
- Feature count: 5039
- Drugs: 295
- Cell lines: 29
- Invalid SMILES pairs kept: 1294

### drug_group_4fold
- weighted_top3_ensemble: Spearman 0.5446, RMSE 2.1413
- xgboost: Spearman 0.5296, RMSE 2.1591
- residual_mlp: Spearman 0.5014, RMSE 2.3047
- flat_mlp: Spearman 0.4987, RMSE 2.2402
- randomforest: Spearman 0.4976, RMSE 2.2164

### random_4fold
- lightgbm: Spearman 0.8396, RMSE 1.2098
- weighted_top3_ensemble: Spearman 0.8281, RMSE 1.2482
- xgboost: Spearman 0.8280, RMSE 1.2562
- randomforest: Spearman 0.7928, RMSE 1.3601
- extratrees: Spearman 0.7805, RMSE 1.4137

### scaffold_group_4fold
- weighted_top3_ensemble: Spearman 0.5578, RMSE 2.1528
- flat_mlp: Spearman 0.5205, RMSE 2.2115
- xgboost: Spearman 0.5139, RMSE 2.2337
- residual_mlp: Spearman 0.5084, RMSE 2.2882
- randomforest: Spearman 0.5026, RMSE 2.2409

## colonstyle_compact_baseline

- Input shape: [7513, 6524]
- Feature count: 4684
- Drugs: 295
- Cell lines: 29
- Invalid SMILES pairs kept: 1294

### drug_group_4fold
- weighted_top3_ensemble: Spearman 0.5418, RMSE 2.1909
- randomforest: Spearman 0.5283, RMSE 2.2065
- extratrees: Spearman 0.5050, RMSE 2.2432
- xgboost: Spearman 0.4893, RMSE 2.2799
- lightgbm: Spearman 0.4832, RMSE 2.3251

### random_4fold
- lightgbm: Spearman 0.8320, RMSE 1.2230
- weighted_top3_ensemble: Spearman 0.8255, RMSE 1.2589
- xgboost: Spearman 0.8225, RMSE 1.2708
- randomforest: Spearman 0.7969, RMSE 1.3643
- extratrees: Spearman 0.7839, RMSE 1.4184

### scaffold_group_4fold
- weighted_top3_ensemble: Spearman 0.5072, RMSE 2.2842
- randomforest: Spearman 0.5028, RMSE 2.2908
- xgboost: Spearman 0.4609, RMSE 2.3582
- extratrees: Spearman 0.4581, RMSE 2.3576
- lightgbm: Spearman 0.4287, RMSE 2.4677

## Variant Policy Comparison

| Split | Best valid-only | Best all-drugs-zero | Delta valid-zero |
|---|---:|---:|---:|
| drug_group_4fold | 0.5823 | 0.5446 | +0.0377 |
| random_4fold | 0.8893 | 0.8396 | +0.0497 |
| scaffold_group_4fold | 0.5775 | 0.5578 | +0.0197 |

## Recommended Input Policy

- Recommended variant: `legacy_rich_valid_smiles_only`
- Basis: highest mean Spearman across drug_group_4fold and scaffold_group_4fold

| Variant | Mean drug/scaffold Spearman |
|---|---:|
| legacy_rich_valid_smiles_only | 0.5799 |
| legacy_rich_all_drugs_zero_smiles | 0.5512 |
| colonstyle_compact_baseline | 0.5245 |

## Extended Metrics For Split Winners

| Split | Variant | Model | Spearman | Pearson | Kendall | RMSE | MAE | R2 | AUROC | AUPRC |
|---|---|---|---:|---:|---:|---:|---:|---:|---:|---:|
| drug_group_4fold | legacy_rich_valid_smiles_only | weighted_top3_ensemble | 0.5823 | 0.6211 | 0.4188 | 2.1716 | 1.6378 | 0.3852 | 0.7996 | 0.6929 |
| random_4fold | legacy_rich_valid_smiles_only | lightgbm | 0.8893 | 0.9252 | 0.7227 | 1.0549 | 0.8105 | 0.8549 | 0.9538 | 0.9163 |
| scaffold_group_4fold | legacy_rich_valid_smiles_only | weighted_top3_ensemble | 0.5775 | 0.5864 | 0.4120 | 2.2475 | 1.6819 | 0.3414 | 0.8010 | 0.6565 |

Binary metrics use `label_binary=1` as sensitive and score `-predicted LN_IC50`.
