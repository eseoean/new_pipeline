# HNSC Hybrid Pipeline Report

- Cancer: Head and neck squamous cell carcinoma (두경부 편평세포암)
- Run ID: `260423_HNSC_V2_raw_lincs_cmap2020_bicr6`
- Raw S3: `s3://say2-4team/HNSC_raw`
- Upload S3: `s3://say2-4team/20260409_eseo/260423_HNSC_V2/raw_lincs_cmap2020_bicr6`
- Fold policy: 4-fold random/drug/scaffold

## Input QC

- Labels: 9358 rows, 39 cell lines, 295 drugs
- Valid SMILES drugs: 186 / 295
- DepMap mapped cell lines: 39 / 39
- CRISPR features: 3500
- LINCS mapped drugs: 38

## LINCS Policy

- Strategy: `cancer_specific_unavailable`
- Main-run eligible: no
- Disease-overlap LINCS cells: none
- Representative cells declared: none
- Matched drugs: 38
- Matched signatures: 114
- Reason: No overlapping disease LINCS cell line was available, so a cancer-specific LINCS main run is unavailable.

## legacy_rich_valid_smiles_only

- Input shape: [6259, 6539]
- Feature count: 5052
- Drugs: 186
- Cell lines: 39
- Invalid SMILES pairs kept: 0

### drug_group_4fold
- weighted_top3_ensemble: Spearman 0.5720, RMSE 2.2716
- randomforest: Spearman 0.5692, RMSE 2.2819
- extratrees: Spearman 0.5528, RMSE 2.2791
- xgboost: Spearman 0.5011, RMSE 2.3972
- flat_mlp: Spearman 0.4915, RMSE 2.3815

### random_4fold
- lightgbm: Spearman 0.9056, RMSE 1.0497
- weighted_top3_ensemble: Spearman 0.9049, RMSE 1.0651
- xgboost: Spearman 0.8969, RMSE 1.1126
- flat_mlp: Spearman 0.8927, RMSE 1.1421
- randomforest: Spearman 0.8921, RMSE 1.1589

### scaffold_group_4fold
- weighted_top3_ensemble: Spearman 0.5765, RMSE 2.4334
- randomforest: Spearman 0.5596, RMSE 2.4636
- extratrees: Spearman 0.5484, RMSE 2.4586
- xgboost: Spearman 0.5211, RMSE 2.5143
- lightgbm: Spearman 0.4681, RMSE 2.6584

## legacy_rich_all_drugs_zero_smiles

- Input shape: [9358, 6539]
- Feature count: 5055
- Drugs: 295
- Cell lines: 39
- Invalid SMILES pairs kept: 3099

### drug_group_4fold
- weighted_top3_ensemble: Spearman 0.5253, RMSE 2.1264
- randomforest: Spearman 0.5175, RMSE 2.1346
- extratrees: Spearman 0.5036, RMSE 2.1591
- xgboost: Spearman 0.4843, RMSE 2.2365
- flat_mlp: Spearman 0.4699, RMSE 2.1797

### random_4fold
- lightgbm: Spearman 0.8207, RMSE 1.2867
- weighted_top3_ensemble: Spearman 0.8136, RMSE 1.3300
- xgboost: Spearman 0.8015, RMSE 1.3691
- flat_mlp: Spearman 0.7764, RMSE 1.4371
- residual_mlp: Spearman 0.7683, RMSE 1.4820

### scaffold_group_4fold
- weighted_top3_ensemble: Spearman 0.5069, RMSE 2.1627
- randomforest: Spearman 0.4936, RMSE 2.2186
- extratrees: Spearman 0.4917, RMSE 2.2137
- residual_mlp: Spearman 0.4620, RMSE 2.2638
- xgboost: Spearman 0.4524, RMSE 2.3471

## colonstyle_compact_baseline

- Input shape: [9358, 6539]
- Feature count: 4700
- Drugs: 295
- Cell lines: 39
- Invalid SMILES pairs kept: 3099

### drug_group_4fold
- randomforest: Spearman 0.5064, RMSE 2.2243
- weighted_top3_ensemble: Spearman 0.4928, RMSE 2.2274
- extratrees: Spearman 0.4711, RMSE 2.2558
- lightgbm: Spearman 0.4404, RMSE 2.3564
- xgboost: Spearman 0.4334, RMSE 2.3536

### random_4fold
- lightgbm: Spearman 0.7930, RMSE 1.3876
- weighted_top3_ensemble: Spearman 0.7863, RMSE 1.4075
- xgboost: Spearman 0.7739, RMSE 1.4514
- residual_mlp: Spearman 0.7571, RMSE 1.4883
- flat_mlp: Spearman 0.7482, RMSE 1.4881

### scaffold_group_4fold
- extratrees: Spearman 0.4653, RMSE 2.3282
- weighted_top3_ensemble: Spearman 0.4644, RMSE 2.3125
- randomforest: Spearman 0.4631, RMSE 2.3399
- lightgbm: Spearman 0.4177, RMSE 2.4411
- xgboost: Spearman 0.4167, RMSE 2.4367

## Variant Policy Comparison

| Split | Best valid-only | Best all-drugs-zero | Delta valid-zero |
|---|---:|---:|---:|
| drug_group_4fold | 0.5720 | 0.5253 | +0.0467 |
| random_4fold | 0.9056 | 0.8207 | +0.0849 |
| scaffold_group_4fold | 0.5765 | 0.5069 | +0.0696 |

## Recommended Input Policy

- Recommended variant: `legacy_rich_valid_smiles_only`
- Basis: highest mean Spearman across drug_group_4fold and scaffold_group_4fold

| Variant | Mean drug/scaffold Spearman |
|---|---:|
| legacy_rich_valid_smiles_only | 0.5743 |
| legacy_rich_all_drugs_zero_smiles | 0.5161 |
| colonstyle_compact_baseline | 0.4859 |

## Extended Metrics For Split Winners

| Split | Variant | Model | Spearman | Pearson | Kendall | RMSE | MAE | R2 | AUROC | AUPRC |
|---|---|---|---:|---:|---:|---:|---:|---:|---:|---:|
| drug_group_4fold | legacy_rich_valid_smiles_only | weighted_top3_ensemble | 0.5720 | 0.6193 | 0.4065 | 2.2716 | 1.7435 | 0.3817 | 0.7951 | 0.7131 |
| random_4fold | legacy_rich_valid_smiles_only | lightgbm | 0.9056 | 0.9322 | 0.7428 | 1.0497 | 0.7832 | 0.8680 | 0.9580 | 0.9315 |
| scaffold_group_4fold | legacy_rich_valid_smiles_only | weighted_top3_ensemble | 0.5765 | 0.5604 | 0.4034 | 2.4334 | 1.8151 | 0.2905 | 0.7863 | 0.6885 |

Binary metrics use `label_binary=1` as sensitive and score `-predicted LN_IC50`.
