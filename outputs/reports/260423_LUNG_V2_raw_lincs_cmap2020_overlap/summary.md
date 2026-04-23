# LUNG Hybrid Pipeline Report

- Cancer: Non-small cell lung cancer (비소세포폐암)
- Run ID: `260423_LUNG_V2_raw_lincs_cmap2020_overlap`
- Raw S3: `s3://say2-4team/Lung_raw`
- Upload S3: `s3://say2-4team/20260409_eseo/260423_LUNG_V2/raw_lincs_cmap2020_overlap`
- Fold policy: 4-fold random/drug/scaffold

## Input QC

- Labels: 19795 rows, 76 cell lines, 295 drugs
- Valid SMILES drugs: 186 / 295
- DepMap mapped cell lines: 77 / 77
- CRISPR features: 3500
- LINCS mapped drugs: 152

## LINCS Policy

- Strategy: `cancer_specific_main`
- Main-run eligible: yes
- Disease-overlap LINCS cells: A549, HCC15, HCC827, NCIH1437, NCIH1563, NCIH1573, NCIH1781, NCIH1975, NCIH838, SKLU1
- Representative cells declared: none
- Matched drugs: 152
- Matched signatures: 5649
- Reason: 10 overlapping disease LINCS cell lines were available (A549, HCC15, HCC827, NCIH1437, NCIH1563, NCIH1573, NCIH1781, NCIH1975, NCIH838, SKLU1), so cancer-specific LINCS can be used as the main setting.

## legacy_rich_valid_smiles_only

- Input shape: [13048, 6574]
- Feature count: 5087
- Drugs: 186
- Cell lines: 76
- Invalid SMILES pairs kept: 0

### drug_group_4fold
- extratrees: Spearman 0.7026, RMSE 2.0512
- weighted_top3_ensemble: Spearman 0.7016, RMSE 2.0279
- randomforest: Spearman 0.6847, RMSE 2.0644
- xgboost: Spearman 0.6612, RMSE 2.1086
- flat_mlp: Spearman 0.6590, RMSE 2.1012

### random_4fold
- weighted_top3_ensemble: Spearman 0.8976, RMSE 1.1575
- lightgbm: Spearman 0.8971, RMSE 1.1508
- flat_mlp: Spearman 0.8907, RMSE 1.1950
- xgboost: Spearman 0.8870, RMSE 1.2132
- residual_mlp: Spearman 0.8765, RMSE 1.2733

### scaffold_group_4fold
- weighted_top3_ensemble: Spearman 0.7161, RMSE 2.0147
- extratrees: Spearman 0.7091, RMSE 2.0232
- randomforest: Spearman 0.7018, RMSE 2.0553
- xgboost: Spearman 0.6843, RMSE 2.0873
- lightgbm: Spearman 0.6693, RMSE 2.1103

## legacy_rich_all_drugs_zero_smiles

- Input shape: [20353, 6574]
- Feature count: 5090
- Drugs: 295
- Cell lines: 76
- Invalid SMILES pairs kept: 7305

### drug_group_4fold
- weighted_top3_ensemble: Spearman 0.6111, RMSE 2.0125
- lightgbm: Spearman 0.6045, RMSE 2.0135
- extratrees: Spearman 0.5876, RMSE 2.0840
- randomforest: Spearman 0.5813, RMSE 2.0831
- xgboost: Spearman 0.5755, RMSE 2.0466

### random_4fold
- lightgbm: Spearman 0.8155, RMSE 1.3599
- weighted_top3_ensemble: Spearman 0.8067, RMSE 1.4108
- xgboost: Spearman 0.7837, RMSE 1.4850
- flat_mlp: Spearman 0.7804, RMSE 1.4918
- residual_mlp: Spearman 0.7779, RMSE 1.5109

### scaffold_group_4fold
- weighted_top3_ensemble: Spearman 0.6030, RMSE 2.0735
- extratrees: Spearman 0.5920, RMSE 2.1005
- randomforest: Spearman 0.5903, RMSE 2.0985
- lightgbm: Spearman 0.5658, RMSE 2.1535
- xgboost: Spearman 0.5626, RMSE 2.1293

## colonstyle_compact_baseline

- Input shape: [20353, 6574]
- Feature count: 4735
- Drugs: 295
- Cell lines: 76
- Invalid SMILES pairs kept: 7305

### drug_group_4fold
- weighted_top3_ensemble: Spearman 0.5635, RMSE 2.1673
- xgboost: Spearman 0.5454, RMSE 2.1949
- lightgbm: Spearman 0.5434, RMSE 2.2063
- randomforest: Spearman 0.5175, RMSE 2.2859
- extratrees: Spearman 0.5088, RMSE 2.2882

### random_4fold
- weighted_top3_ensemble: Spearman 0.7857, RMSE 1.4775
- lightgbm: Spearman 0.7812, RMSE 1.4939
- flat_mlp: Spearman 0.7758, RMSE 1.5054
- residual_mlp: Spearman 0.7726, RMSE 1.5296
- xgboost: Spearman 0.7605, RMSE 1.5635

### scaffold_group_4fold
- weighted_top3_ensemble: Spearman 0.5522, RMSE 2.1876
- lightgbm: Spearman 0.5387, RMSE 2.2252
- xgboost: Spearman 0.5347, RMSE 2.2133
- randomforest: Spearman 0.4988, RMSE 2.3248
- extratrees: Spearman 0.4868, RMSE 2.3559

## Variant Policy Comparison

| Split | Best valid-only | Best all-drugs-zero | Delta valid-zero |
|---|---:|---:|---:|
| drug_group_4fold | 0.7026 | 0.6111 | +0.0915 |
| random_4fold | 0.8976 | 0.8155 | +0.0821 |
| scaffold_group_4fold | 0.7161 | 0.6030 | +0.1131 |

## Recommended Input Policy

- Recommended variant: `legacy_rich_valid_smiles_only`
- Basis: highest mean Spearman across drug_group_4fold and scaffold_group_4fold

| Variant | Mean drug/scaffold Spearman |
|---|---:|
| legacy_rich_valid_smiles_only | 0.7094 |
| legacy_rich_all_drugs_zero_smiles | 0.6071 |
| colonstyle_compact_baseline | 0.5578 |

## Extended Metrics For Split Winners

| Split | Variant | Model | Spearman | Pearson | Kendall | RMSE | MAE | R2 | AUROC | AUPRC |
|---|---|---|---:|---:|---:|---:|---:|---:|---:|---:|
| drug_group_4fold | legacy_rich_valid_smiles_only | extratrees | 0.7026 | 0.7233 | 0.5198 | 2.0512 | 1.5502 | 0.5167 | 0.8602 | 0.8070 |
| random_4fold | legacy_rich_valid_smiles_only | weighted_top3_ensemble | 0.8976 | 0.9214 | 0.7305 | 1.1575 | 0.8796 | 0.8461 | 0.9509 | 0.9208 |
| scaffold_group_4fold | legacy_rich_valid_smiles_only | weighted_top3_ensemble | 0.7161 | 0.7457 | 0.5323 | 2.0147 | 1.5186 | 0.5338 | 0.8643 | 0.8103 |

Binary metrics use `label_binary=1` as sensitive and score `-predicted LN_IC50`.
