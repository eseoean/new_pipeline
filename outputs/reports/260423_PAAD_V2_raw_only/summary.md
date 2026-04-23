# PAAD Hybrid Pipeline Report

- Cancer: Pancreatic adenocarcinoma (췌장암)
- Run ID: `260423_PAAD_V2_raw_only`
- Raw S3: `s3://say2-4team/PAAD_raw`
- Upload S3: `s3://say2-4team/20260409_eseo/260423_PAAD_V2/raw_only`
- Fold policy: 4-fold random/drug/scaffold

## Input QC

- Labels: 7513 rows, 29 cell lines, 295 drugs
- Valid SMILES drugs: 160 / 295
- DepMap mapped cell lines: 29 / 29
- CRISPR features: 3500
- LINCS mapped drugs: 120

## legacy_rich_valid_smiles_only

- Input shape: [4110, 6524]
- Feature count: 5036
- Drugs: 160
- Cell lines: 29
- Invalid SMILES pairs kept: 0

### drug_group_4fold
- weighted_top3_ensemble: Spearman 0.7123, RMSE 2.0122
- extratrees: Spearman 0.7060, RMSE 2.0127
- randomforest: Spearman 0.7027, RMSE 2.0248
- xgboost: Spearman 0.6739, RMSE 2.1145
- lightgbm: Spearman 0.6431, RMSE 2.2034

### random_4fold
- lightgbm: Spearman 0.9127, RMSE 1.0665
- weighted_top3_ensemble: Spearman 0.9118, RMSE 1.0802
- xgboost: Spearman 0.9080, RMSE 1.0984
- randomforest: Spearman 0.9049, RMSE 1.1353
- extratrees: Spearman 0.8979, RMSE 1.1942

### scaffold_group_4fold
- weighted_top3_ensemble: Spearman 0.7061, RMSE 2.0374
- lightgbm: Spearman 0.6967, RMSE 2.0697
- xgboost: Spearman 0.6893, RMSE 2.0913
- extratrees: Spearman 0.6756, RMSE 2.1180
- randomforest: Spearman 0.6721, RMSE 2.1406

## legacy_rich_all_drugs_zero_smiles

- Input shape: [7513, 6524]
- Feature count: 5039
- Drugs: 295
- Cell lines: 29
- Invalid SMILES pairs kept: 3403

### drug_group_4fold
- weighted_top3_ensemble: Spearman 0.6058, RMSE 1.9565
- xgboost: Spearman 0.5955, RMSE 1.9808
- lightgbm: Spearman 0.5875, RMSE 2.0248
- randomforest: Spearman 0.5836, RMSE 1.9720
- extratrees: Spearman 0.5671, RMSE 2.0117

### random_4fold
- lightgbm: Spearman 0.7988, RMSE 1.3350
- xgboost: Spearman 0.7826, RMSE 1.3978
- weighted_top3_ensemble: Spearman 0.7823, RMSE 1.3895
- randomforest: Spearman 0.7226, RMSE 1.5402
- flat_mlp: Spearman 0.7180, RMSE 1.5785

### scaffold_group_4fold
- weighted_top3_ensemble: Spearman 0.5631, RMSE 2.0223
- randomforest: Spearman 0.5541, RMSE 2.0385
- xgboost: Spearman 0.5486, RMSE 2.0795
- extratrees: Spearman 0.5358, RMSE 2.0797
- residual_mlp: Spearman 0.5349, RMSE 2.1333

## colonstyle_compact_baseline

- Input shape: [7513, 6524]
- Feature count: 4684
- Drugs: 295
- Cell lines: 29
- Invalid SMILES pairs kept: 3403

### drug_group_4fold
- weighted_top3_ensemble: Spearman 0.5422, RMSE 2.0794
- xgboost: Spearman 0.5310, RMSE 2.0940
- lightgbm: Spearman 0.5252, RMSE 2.1108
- randomforest: Spearman 0.5110, RMSE 2.2000
- extratrees: Spearman 0.5103, RMSE 2.2027

### random_4fold
- lightgbm: Spearman 0.7539, RMSE 1.4743
- xgboost: Spearman 0.7436, RMSE 1.5060
- weighted_top3_ensemble: Spearman 0.7417, RMSE 1.5061
- randomforest: Spearman 0.7023, RMSE 1.6132
- flat_mlp: Spearman 0.6983, RMSE 1.6189

### scaffold_group_4fold
- weighted_top3_ensemble: Spearman 0.5351, RMSE 2.0994
- xgboost: Spearman 0.5199, RMSE 2.1030
- lightgbm: Spearman 0.5197, RMSE 2.1359
- randomforest: Spearman 0.4998, RMSE 2.2431
- extratrees: Spearman 0.4808, RMSE 2.2868

## Variant Policy Comparison

| Split | Best valid-only | Best all-drugs-zero | Delta valid-zero |
|---|---:|---:|---:|
| drug_group_4fold | 0.7123 | 0.6058 | +0.1065 |
| random_4fold | 0.9127 | 0.7988 | +0.1139 |
| scaffold_group_4fold | 0.7061 | 0.5631 | +0.1430 |

## Recommended Input Policy

- Recommended variant: `legacy_rich_valid_smiles_only`
- Basis: highest mean Spearman across drug_group_4fold and scaffold_group_4fold

| Variant | Mean drug/scaffold Spearman |
|---|---:|
| legacy_rich_valid_smiles_only | 0.7092 |
| legacy_rich_all_drugs_zero_smiles | 0.5844 |
| colonstyle_compact_baseline | 0.5387 |

## Extended Metrics For Split Winners

| Split | Variant | Model | Spearman | Pearson | Kendall | RMSE | MAE | R2 | AUROC | AUPRC |
|---|---|---|---:|---:|---:|---:|---:|---:|---:|---:|
| drug_group_4fold | legacy_rich_valid_smiles_only | weighted_top3_ensemble | 0.7123 | 0.7490 | 0.5279 | 2.0122 | 1.5434 | 0.5550 | 0.8602 | 0.8142 |
| random_4fold | legacy_rich_valid_smiles_only | lightgbm | 0.9127 | 0.9356 | 0.7530 | 1.0665 | 0.8184 | 0.8750 | 0.9623 | 0.9426 |
| scaffold_group_4fold | legacy_rich_valid_smiles_only | weighted_top3_ensemble | 0.7061 | 0.7406 | 0.5247 | 2.0374 | 1.5331 | 0.5438 | 0.8602 | 0.8145 |

Binary metrics use `label_binary=1` as sensitive and score `-predicted LN_IC50`.
