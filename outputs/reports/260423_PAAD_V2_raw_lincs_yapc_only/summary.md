# PAAD Hybrid Pipeline Report

- Cancer: Pancreatic adenocarcinoma (췌장암)
- Run ID: `260423_PAAD_V2_raw_lincs_yapc_only`
- Raw S3: `s3://say2-4team/PAAD_raw`
- Upload S3: `s3://say2-4team/20260409_eseo/260423_PAAD_V2/raw_lincs_yapc_only`
- Fold policy: 4-fold random/drug/scaffold

## Input QC

- Labels: 7513 rows, 29 cell lines, 295 drugs
- Valid SMILES drugs: 160 / 295
- DepMap mapped cell lines: 29 / 29
- CRISPR features: 3500
- LINCS mapped drugs: 85

## legacy_rich_valid_smiles_only

- Input shape: [4110, 6524]
- Feature count: 5036
- Drugs: 160
- Cell lines: 29
- Invalid SMILES pairs kept: 0

### drug_group_4fold
- weighted_top3_ensemble: Spearman 0.6410, RMSE 2.1689
- randomforest: Spearman 0.6332, RMSE 2.1976
- extratrees: Spearman 0.6145, RMSE 2.2330
- residual_mlp: Spearman 0.6056, RMSE 2.2870
- flat_mlp: Spearman 0.6051, RMSE 2.2047

### random_4fold
- lightgbm: Spearman 0.9110, RMSE 1.0817
- weighted_top3_ensemble: Spearman 0.9103, RMSE 1.0930
- xgboost: Spearman 0.9076, RMSE 1.1039
- randomforest: Spearman 0.9014, RMSE 1.1559
- extratrees: Spearman 0.8978, RMSE 1.1974

### scaffold_group_4fold
- randomforest: Spearman 0.6298, RMSE 2.2593
- extratrees: Spearman 0.6215, RMSE 2.2508
- weighted_top3_ensemble: Spearman 0.6193, RMSE 2.2708
- xgboost: Spearman 0.5506, RMSE 2.4421
- residual_mlp: Spearman 0.5493, RMSE 2.6228

## legacy_rich_all_drugs_zero_smiles

- Input shape: [7513, 6524]
- Feature count: 5039
- Drugs: 295
- Cell lines: 29
- Invalid SMILES pairs kept: 3403

### drug_group_4fold
- weighted_top3_ensemble: Spearman 0.5536, RMSE 2.0374
- randomforest: Spearman 0.5419, RMSE 2.0600
- xgboost: Spearman 0.5417, RMSE 2.0567
- lightgbm: Spearman 0.5328, RMSE 2.1241
- extratrees: Spearman 0.5263, RMSE 2.1076

### random_4fold
- lightgbm: Spearman 0.7973, RMSE 1.3384
- weighted_top3_ensemble: Spearman 0.7861, RMSE 1.3920
- xgboost: Spearman 0.7821, RMSE 1.3971
- flat_mlp: Spearman 0.7170, RMSE 1.5887
- randomforest: Spearman 0.7167, RMSE 1.5543

### scaffold_group_4fold
- weighted_top3_ensemble: Spearman 0.5443, RMSE 2.1057
- randomforest: Spearman 0.5196, RMSE 2.1342
- residual_mlp: Spearman 0.5195, RMSE 2.2046
- flat_mlp: Spearman 0.4971, RMSE 2.1747
- extratrees: Spearman 0.4801, RMSE 2.2039

## colonstyle_compact_baseline

- Input shape: [7513, 6524]
- Feature count: 4684
- Drugs: 295
- Cell lines: 29
- Invalid SMILES pairs kept: 3403

### drug_group_4fold
- weighted_top3_ensemble: Spearman 0.5369, RMSE 2.1474
- xgboost: Spearman 0.5183, RMSE 2.1298
- randomforest: Spearman 0.5118, RMSE 2.2210
- extratrees: Spearman 0.5057, RMSE 2.2276
- lightgbm: Spearman 0.5030, RMSE 2.1727

### random_4fold
- lightgbm: Spearman 0.7547, RMSE 1.4738
- weighted_top3_ensemble: Spearman 0.7479, RMSE 1.4902
- xgboost: Spearman 0.7437, RMSE 1.5057
- flat_mlp: Spearman 0.7032, RMSE 1.5990
- randomforest: Spearman 0.7021, RMSE 1.6110

### scaffold_group_4fold
- weighted_top3_ensemble: Spearman 0.5199, RMSE 2.1782
- lightgbm: Spearman 0.5086, RMSE 2.1763
- xgboost: Spearman 0.4997, RMSE 2.2019
- randomforest: Spearman 0.4873, RMSE 2.2967
- extratrees: Spearman 0.4718, RMSE 2.3283

## Variant Policy Comparison

| Split | Best valid-only | Best all-drugs-zero | Delta valid-zero |
|---|---:|---:|---:|
| drug_group_4fold | 0.6410 | 0.5536 | +0.0874 |
| random_4fold | 0.9110 | 0.7973 | +0.1137 |
| scaffold_group_4fold | 0.6298 | 0.5443 | +0.0855 |

## Recommended Input Policy

- Recommended variant: `legacy_rich_valid_smiles_only`
- Basis: highest mean Spearman across drug_group_4fold and scaffold_group_4fold

| Variant | Mean drug/scaffold Spearman |
|---|---:|
| legacy_rich_valid_smiles_only | 0.6354 |
| legacy_rich_all_drugs_zero_smiles | 0.5489 |
| colonstyle_compact_baseline | 0.5284 |

## Extended Metrics For Split Winners

| Split | Variant | Model | Spearman | Pearson | Kendall | RMSE | MAE | R2 | AUROC | AUPRC |
|---|---|---|---:|---:|---:|---:|---:|---:|---:|---:|
| drug_group_4fold | legacy_rich_valid_smiles_only | weighted_top3_ensemble | 0.6410 | 0.6955 | 0.4701 | 2.1689 | 1.6504 | 0.4830 | 0.8254 | 0.7873 |
| random_4fold | legacy_rich_valid_smiles_only | lightgbm | 0.9110 | 0.9337 | 0.7501 | 1.0817 | 0.8307 | 0.8714 | 0.9617 | 0.9415 |
| scaffold_group_4fold | legacy_rich_valid_smiles_only | randomforest | 0.6298 | 0.6708 | 0.4595 | 2.2593 | 1.7052 | 0.4390 | 0.8225 | 0.7858 |

Binary metrics use `label_binary=1` as sensitive and score `-predicted LN_IC50`.
