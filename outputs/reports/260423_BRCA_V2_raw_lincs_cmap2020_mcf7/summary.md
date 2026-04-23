# BRCA Hybrid Pipeline Report

- Cancer: Breast invasive carcinoma (유방암)
- Run ID: `260423_BRCA_V2_raw_lincs_cmap2020_mcf7`
- Raw S3: `s3://say2-4team/oringinal_raw`
- Upload S3: `s3://say2-4team/20260409_eseo/260423_BRCA_V2/raw_lincs_cmap2020_mcf7`
- Fold policy: 4-fold random/drug/scaffold

## Input QC

- Labels: 13106 rows, 51 cell lines, 295 drugs
- Valid SMILES drugs: 186 / 295
- DepMap mapped cell lines: 51 / 51
- CRISPR features: 3500
- LINCS mapped drugs: 148

## LINCS Policy

- Strategy: `cancer_specific_main`
- Main-run eligible: yes
- Disease-overlap LINCS cells: BT20, BT474, HS578T, MCF7, MDAMB231, MDAMB468, T47D
- Representative cells declared: MCF7
- Matched drugs: 148
- Matched signatures: 6224
- Reason: 7 overlapping disease LINCS cell lines were available (BT20, BT474, HS578T, MCF7, MDAMB231, MDAMB468, T47D), so cancer-specific LINCS can be used as the main setting.

## legacy_rich_valid_smiles_only

- Input shape: [8409, 6551]
- Feature count: 5064
- Drugs: 186
- Cell lines: 51
- Invalid SMILES pairs kept: 0

### drug_group_4fold
- weighted_top3_ensemble: Spearman 0.6708, RMSE 2.1150
- extratrees: Spearman 0.6628, RMSE 2.0844
- randomforest: Spearman 0.6472, RMSE 2.1337
- xgboost: Spearman 0.6322, RMSE 2.2535
- lightgbm: Spearman 0.6255, RMSE 2.2983

### random_4fold
- weighted_top3_ensemble: Spearman 0.8699, RMSE 1.2388
- lightgbm: Spearman 0.8680, RMSE 1.2361
- xgboost: Spearman 0.8584, RMSE 1.2857
- flat_mlp: Spearman 0.8576, RMSE 1.2940
- randomforest: Spearman 0.8284, RMSE 1.4125

### scaffold_group_4fold
- weighted_top3_ensemble: Spearman 0.6821, RMSE 2.0001
- extratrees: Spearman 0.6619, RMSE 2.0514
- flat_mlp: Spearman 0.6451, RMSE 2.0507
- xgboost: Spearman 0.6306, RMSE 2.0873
- randomforest: Spearman 0.6274, RMSE 2.1383

## legacy_rich_all_drugs_zero_smiles

- Input shape: [13106, 6551]
- Feature count: 5067
- Drugs: 295
- Cell lines: 51
- Invalid SMILES pairs kept: 4697

### drug_group_4fold
- weighted_top3_ensemble: Spearman 0.5724, RMSE 2.0626
- xgboost: Spearman 0.5621, RMSE 2.1101
- lightgbm: Spearman 0.5541, RMSE 2.1272
- extratrees: Spearman 0.5407, RMSE 2.0900
- randomforest: Spearman 0.5405, RMSE 2.1048

### random_4fold
- lightgbm: Spearman 0.7830, RMSE 1.4211
- weighted_top3_ensemble: Spearman 0.7779, RMSE 1.4583
- xgboost: Spearman 0.7603, RMSE 1.5093
- flat_mlp: Spearman 0.7495, RMSE 1.5470
- residual_mlp: Spearman 0.7358, RMSE 1.5970

### scaffold_group_4fold
- weighted_top3_ensemble: Spearman 0.5495, RMSE 2.0760
- lightgbm: Spearman 0.5329, RMSE 2.1256
- xgboost: Spearman 0.5252, RMSE 2.1401
- extratrees: Spearman 0.5238, RMSE 2.1149
- randomforest: Spearman 0.5233, RMSE 2.1099

## colonstyle_compact_baseline

- Input shape: [13106, 6551]
- Feature count: 4712
- Drugs: 295
- Cell lines: 51
- Invalid SMILES pairs kept: 4697

### drug_group_4fold
- weighted_top3_ensemble: Spearman 0.5243, RMSE 2.1526
- xgboost: Spearman 0.4996, RMSE 2.1971
- randomforest: Spearman 0.4980, RMSE 2.1998
- extratrees: Spearman 0.4885, RMSE 2.2122
- lightgbm: Spearman 0.4862, RMSE 2.2302

### random_4fold
- weighted_top3_ensemble: Spearman 0.7507, RMSE 1.5385
- lightgbm: Spearman 0.7505, RMSE 1.5382
- xgboost: Spearman 0.7341, RMSE 1.5874
- flat_mlp: Spearman 0.7302, RMSE 1.5726
- residual_mlp: Spearman 0.7287, RMSE 1.6016

### scaffold_group_4fold
- weighted_top3_ensemble: Spearman 0.5133, RMSE 2.1773
- xgboost: Spearman 0.5069, RMSE 2.2278
- lightgbm: Spearman 0.4955, RMSE 2.2011
- randomforest: Spearman 0.4525, RMSE 2.2743
- extratrees: Spearman 0.4430, RMSE 2.2912

## Variant Policy Comparison

| Split | Best valid-only | Best all-drugs-zero | Delta valid-zero |
|---|---:|---:|---:|
| drug_group_4fold | 0.6708 | 0.5724 | +0.0985 |
| random_4fold | 0.8699 | 0.7830 | +0.0869 |
| scaffold_group_4fold | 0.6821 | 0.5495 | +0.1326 |

## Recommended Input Policy

- Recommended variant: `legacy_rich_valid_smiles_only`
- Basis: highest mean Spearman across drug_group_4fold and scaffold_group_4fold

| Variant | Mean drug/scaffold Spearman |
|---|---:|
| legacy_rich_valid_smiles_only | 0.6764 |
| legacy_rich_all_drugs_zero_smiles | 0.5609 |
| colonstyle_compact_baseline | 0.5188 |

## Extended Metrics For Split Winners

| Split | Variant | Model | Spearman | Pearson | Kendall | RMSE | MAE | R2 | AUROC | AUPRC |
|---|---|---|---:|---:|---:|---:|---:|---:|---:|---:|
| drug_group_4fold | legacy_rich_valid_smiles_only | weighted_top3_ensemble | 0.6708 | 0.6766 | 0.4910 | 2.1150 | 1.5818 | 0.4527 | 0.8426 | 0.7653 |
| random_4fold | legacy_rich_valid_smiles_only | weighted_top3_ensemble | 0.8699 | 0.9031 | 0.6957 | 1.2388 | 0.9246 | 0.8123 | 0.9401 | 0.9006 |
| scaffold_group_4fold | legacy_rich_valid_smiles_only | weighted_top3_ensemble | 0.6821 | 0.7267 | 0.5036 | 2.0001 | 1.5262 | 0.5106 | 0.8437 | 0.7827 |

Binary metrics use `label_binary=1` as sensitive and score `-predicted LN_IC50`.
