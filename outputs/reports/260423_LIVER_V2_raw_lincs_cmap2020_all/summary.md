# LIHC Hybrid Pipeline Report

- Cancer: Liver hepatocellular carcinoma (간암)
- Run ID: `260423_LIVER_V2_raw_lincs_cmap2020_all`
- Raw S3: `s3://say2-4team/Liver_raw`
- Upload S3: `s3://say2-4team/20260409_eseo/260423_LIVER_V2/raw_lincs_cmap2020_all`
- Fold policy: 4-fold random/drug/scaffold

## Input QC

- Labels: 4164 rows, 15 cell lines, 295 drugs
- Valid SMILES drugs: 186 / 295
- DepMap mapped cell lines: 15 / 15
- CRISPR features: 3500
- LINCS mapped drugs: 176

## LINCS Policy

- Strategy: `all_cell_main`
- Main-run eligible: no
- Disease-overlap LINCS cells: HUH7, JHH7
- Representative cells declared: none
- Matched drugs: 176
- Matched signatures: 81207
- Reason: This run keeps all-cell LINCS although 2 overlapping disease LINCS cell lines were available (HUH7, JHH7). In that situation, a disease-restricted LINCS main run should usually be preferred.

## legacy_rich_valid_smiles_only

- Input shape: [2590, 6511]
- Feature count: 5024
- Drugs: 186
- Cell lines: 15
- Invalid SMILES pairs kept: 0

### drug_group_4fold
- extratrees: Spearman 0.6910, RMSE 1.8775
- weighted_top3_ensemble: Spearman 0.6878, RMSE 1.8799
- randomforest: Spearman 0.6776, RMSE 1.9068
- xgboost: Spearman 0.6616, RMSE 1.9675
- lightgbm: Spearman 0.6497, RMSE 2.0170

### random_4fold
- lightgbm: Spearman 0.8752, RMSE 1.1573
- weighted_top3_ensemble: Spearman 0.8725, RMSE 1.1700
- xgboost: Spearman 0.8715, RMSE 1.1674
- flat_mlp: Spearman 0.8391, RMSE 1.3178
- randomforest: Spearman 0.8342, RMSE 1.2994

### scaffold_group_4fold
- weighted_top3_ensemble: Spearman 0.6712, RMSE 2.0228
- extratrees: Spearman 0.6680, RMSE 2.0189
- randomforest: Spearman 0.6567, RMSE 2.0564
- xgboost: Spearman 0.6478, RMSE 2.0700
- residual_mlp: Spearman 0.6435, RMSE 2.0637

## legacy_rich_all_drugs_zero_smiles

- Input shape: [4164, 6511]
- Feature count: 5027
- Drugs: 295
- Cell lines: 15
- Invalid SMILES pairs kept: 1574

### drug_group_4fold
- weighted_top3_ensemble: Spearman 0.5669, RMSE 2.0072
- extratrees: Spearman 0.5563, RMSE 2.0069
- randomforest: Spearman 0.5500, RMSE 2.0232
- xgboost: Spearman 0.5401, RMSE 2.0922
- lightgbm: Spearman 0.5269, RMSE 2.0905

### random_4fold
- lightgbm: Spearman 0.7638, RMSE 1.4180
- xgboost: Spearman 0.7617, RMSE 1.4387
- weighted_top3_ensemble: Spearman 0.7537, RMSE 1.4511
- randomforest: Spearman 0.6906, RMSE 1.6021
- extratrees: Spearman 0.6792, RMSE 1.6456

### scaffold_group_4fold
- weighted_top3_ensemble: Spearman 0.5834, RMSE 1.9378
- xgboost: Spearman 0.5702, RMSE 1.9998
- extratrees: Spearman 0.5671, RMSE 1.9560
- randomforest: Spearman 0.5643, RMSE 1.9750
- lightgbm: Spearman 0.5463, RMSE 2.0577

## colonstyle_compact_baseline

- Input shape: [4164, 6511]
- Feature count: 4672
- Drugs: 295
- Cell lines: 15
- Invalid SMILES pairs kept: 1574

### drug_group_4fold
- xgboost: Spearman 0.4755, RMSE 2.1778
- weighted_top3_ensemble: Spearman 0.4743, RMSE 2.1613
- randomforest: Spearman 0.4451, RMSE 2.2311
- lightgbm: Spearman 0.4403, RMSE 2.2177
- extratrees: Spearman 0.4305, RMSE 2.2572

### random_4fold
- lightgbm: Spearman 0.7314, RMSE 1.5385
- weighted_top3_ensemble: Spearman 0.7281, RMSE 1.5539
- xgboost: Spearman 0.7250, RMSE 1.5553
- flat_mlp: Spearman 0.6550, RMSE 1.6978
- randomforest: Spearman 0.6494, RMSE 1.7113

### scaffold_group_4fold
- weighted_top3_ensemble: Spearman 0.4789, RMSE 2.1464
- xgboost: Spearman 0.4752, RMSE 2.1556
- lightgbm: Spearman 0.4638, RMSE 2.2397
- randomforest: Spearman 0.4281, RMSE 2.2528
- extratrees: Spearman 0.4051, RMSE 2.2802

## Variant Policy Comparison

| Split | Best valid-only | Best all-drugs-zero | Delta valid-zero |
|---|---:|---:|---:|
| drug_group_4fold | 0.6910 | 0.5669 | +0.1241 |
| random_4fold | 0.8752 | 0.7638 | +0.1115 |
| scaffold_group_4fold | 0.6712 | 0.5834 | +0.0878 |

## Recommended Input Policy

- Recommended variant: `legacy_rich_valid_smiles_only`
- Basis: highest mean Spearman across drug_group_4fold and scaffold_group_4fold

| Variant | Mean drug/scaffold Spearman |
|---|---:|
| legacy_rich_valid_smiles_only | 0.6811 |
| legacy_rich_all_drugs_zero_smiles | 0.5752 |
| colonstyle_compact_baseline | 0.4772 |

## Extended Metrics For Split Winners

| Split | Variant | Model | Spearman | Pearson | Kendall | RMSE | MAE | R2 | AUROC | AUPRC |
|---|---|---|---:|---:|---:|---:|---:|---:|---:|---:|
| drug_group_4fold | legacy_rich_valid_smiles_only | extratrees | 0.6910 | 0.7310 | 0.5085 | 1.8775 | 1.4205 | 0.5310 | 0.8567 | 0.7854 |
| random_4fold | legacy_rich_valid_smiles_only | lightgbm | 0.8752 | 0.9065 | 0.7013 | 1.1573 | 0.8761 | 0.8218 | 0.9482 | 0.9073 |
| scaffold_group_4fold | legacy_rich_valid_smiles_only | weighted_top3_ensemble | 0.6712 | 0.6898 | 0.4913 | 2.0228 | 1.5139 | 0.4555 | 0.8523 | 0.7735 |

Binary metrics use `label_binary=1` as sensitive and score `-predicted LN_IC50`.
