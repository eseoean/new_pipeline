# PAAD Hybrid Pipeline Report

- Cancer: Pancreatic adenocarcinoma (췌장암)
- Run ID: `260423_PAAD_V2_original_core`
- Raw S3: `s3://say2-4team/PAAD_raw`
- Upload S3: `s3://say2-4team/20260409_eseo/260423_PAAD_V2/original_core`
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
- weighted_top3_ensemble: Spearman 0.5768, RMSE 2.1862
- flat_mlp: Spearman 0.5560, RMSE 2.2033
- randomforest: Spearman 0.5380, RMSE 2.2632
- xgboost: Spearman 0.5342, RMSE 2.2750
- extratrees: Spearman 0.5235, RMSE 2.2564

### random_4fold
- weighted_top3_ensemble: Spearman 0.8887, RMSE 1.0673
- lightgbm: Spearman 0.8875, RMSE 1.0652
- xgboost: Spearman 0.8810, RMSE 1.1135
- flat_mlp: Spearman 0.8713, RMSE 1.1437
- randomforest: Spearman 0.8713, RMSE 1.1597

### scaffold_group_4fold
- weighted_top3_ensemble: Spearman 0.5758, RMSE 2.2477
- xgboost: Spearman 0.5587, RMSE 2.2801
- lightgbm: Spearman 0.5537, RMSE 2.3139
- randomforest: Spearman 0.5493, RMSE 2.2896
- extratrees: Spearman 0.5361, RMSE 2.3124

## legacy_rich_all_drugs_zero_smiles

- Input shape: [7513, 6524]
- Feature count: 5039
- Drugs: 295
- Cell lines: 29
- Invalid SMILES pairs kept: 1294

### drug_group_4fold
- weighted_top3_ensemble: Spearman 0.5470, RMSE 2.1524
- xgboost: Spearman 0.5220, RMSE 2.1905
- flat_mlp: Spearman 0.5118, RMSE 2.2189
- residual_mlp: Spearman 0.5000, RMSE 2.3213
- randomforest: Spearman 0.4869, RMSE 2.2379

### random_4fold
- lightgbm: Spearman 0.8388, RMSE 1.2083
- weighted_top3_ensemble: Spearman 0.8288, RMSE 1.2484
- xgboost: Spearman 0.8279, RMSE 1.2612
- randomforest: Spearman 0.7950, RMSE 1.3569
- flat_mlp: Spearman 0.7847, RMSE 1.4102

### scaffold_group_4fold
- weighted_top3_ensemble: Spearman 0.5530, RMSE 2.1723
- xgboost: Spearman 0.5229, RMSE 2.2358
- flat_mlp: Spearman 0.5140, RMSE 2.2252
- residual_mlp: Spearman 0.5017, RMSE 2.3162
- randomforest: Spearman 0.4964, RMSE 2.2521

## colonstyle_compact_baseline

- Input shape: [7513, 6524]
- Feature count: 4684
- Drugs: 295
- Cell lines: 29
- Invalid SMILES pairs kept: 1294

### drug_group_4fold
- weighted_top3_ensemble: Spearman 0.5415, RMSE 2.1891
- randomforest: Spearman 0.5383, RMSE 2.1914
- extratrees: Spearman 0.5011, RMSE 2.2501
- xgboost: Spearman 0.4856, RMSE 2.2769
- lightgbm: Spearman 0.4810, RMSE 2.3274

### random_4fold
- lightgbm: Spearman 0.8321, RMSE 1.2228
- weighted_top3_ensemble: Spearman 0.8253, RMSE 1.2627
- xgboost: Spearman 0.8230, RMSE 1.2717
- randomforest: Spearman 0.7950, RMSE 1.3745
- extratrees: Spearman 0.7840, RMSE 1.4168

### scaffold_group_4fold
- weighted_top3_ensemble: Spearman 0.5082, RMSE 2.2755
- randomforest: Spearman 0.4957, RMSE 2.3039
- xgboost: Spearman 0.4663, RMSE 2.3163
- extratrees: Spearman 0.4540, RMSE 2.3634
- lightgbm: Spearman 0.4289, RMSE 2.4669

## Variant Policy Comparison

| Split | Best valid-only | Best all-drugs-zero | Delta valid-zero |
|---|---:|---:|---:|
| drug_group_4fold | 0.5768 | 0.5470 | +0.0298 |
| random_4fold | 0.8887 | 0.8388 | +0.0499 |
| scaffold_group_4fold | 0.5758 | 0.5530 | +0.0228 |

## Recommended Input Policy

- Recommended variant: `legacy_rich_valid_smiles_only`
- Basis: highest mean Spearman across drug_group_4fold and scaffold_group_4fold

| Variant | Mean drug/scaffold Spearman |
|---|---:|
| legacy_rich_valid_smiles_only | 0.5763 |
| legacy_rich_all_drugs_zero_smiles | 0.5500 |
| colonstyle_compact_baseline | 0.5249 |

## Extended Metrics For Split Winners

| Split | Variant | Model | Spearman | Pearson | Kendall | RMSE | MAE | R2 | AUROC | AUPRC |
|---|---|---|---:|---:|---:|---:|---:|---:|---:|---:|
| drug_group_4fold | legacy_rich_valid_smiles_only | weighted_top3_ensemble | 0.5768 | 0.6142 | 0.4145 | 2.1862 | 1.6479 | 0.3769 | 0.7974 | 0.6891 |
| random_4fold | legacy_rich_valid_smiles_only | weighted_top3_ensemble | 0.8887 | 0.9249 | 0.7212 | 1.0673 | 0.8272 | 0.8515 | 0.9550 | 0.9185 |
| scaffold_group_4fold | legacy_rich_valid_smiles_only | weighted_top3_ensemble | 0.5758 | 0.5866 | 0.4104 | 2.2477 | 1.6866 | 0.3413 | 0.8001 | 0.6574 |

Binary metrics use `label_binary=1` as sensitive and score `-predicted LN_IC50`.
