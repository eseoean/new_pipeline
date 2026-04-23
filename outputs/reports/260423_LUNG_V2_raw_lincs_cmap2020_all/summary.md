# LUNG Hybrid Pipeline Report

- Cancer: Non-small cell lung cancer (비소세포폐암)
- Run ID: `260423_LUNG_V2_raw_lincs_cmap2020_all`
- Raw S3: `s3://say2-4team/Lung_raw`
- Upload S3: `s3://say2-4team/20260409_eseo/260423_LUNG_V2/raw_lincs_cmap2020_all`
- Fold policy: 4-fold random/drug/scaffold

## Input QC

- Labels: 19795 rows, 76 cell lines, 295 drugs
- Valid SMILES drugs: 186 / 295
- DepMap mapped cell lines: 77 / 77
- CRISPR features: 3500
- LINCS mapped drugs: 176

## LINCS Policy

- Strategy: `all_cell_main`
- Main-run eligible: no
- Disease-overlap LINCS cells: A549, HCC15, HCC827, NCIH1437, NCIH1563, NCIH1573, NCIH1781, NCIH1975, NCIH838, SKLU1
- Representative cells declared: none
- Matched drugs: 176
- Matched signatures: 81207
- Reason: This run keeps all-cell LINCS although 10 overlapping disease LINCS cell lines were available (A549, HCC15, HCC827, NCIH1437, NCIH1563, NCIH1573, NCIH1781, NCIH1975, NCIH838, SKLU1). In that situation, a disease-restricted LINCS main run should usually be preferred.

## legacy_rich_valid_smiles_only

- Input shape: [13048, 6574]
- Feature count: 5087
- Drugs: 186
- Cell lines: 76
- Invalid SMILES pairs kept: 0

### drug_group_4fold
- extratrees: Spearman 0.7430, RMSE 1.9388
- weighted_top3_ensemble: Spearman 0.7385, RMSE 1.9555
- randomforest: Spearman 0.7252, RMSE 1.9815
- xgboost: Spearman 0.6996, RMSE 2.0736
- flat_mlp: Spearman 0.6917, RMSE 2.0525

### random_4fold
- weighted_top3_ensemble: Spearman 0.8978, RMSE 1.1544
- lightgbm: Spearman 0.8965, RMSE 1.1516
- flat_mlp: Spearman 0.8927, RMSE 1.1890
- xgboost: Spearman 0.8880, RMSE 1.2077
- residual_mlp: Spearman 0.8768, RMSE 1.2755

### scaffold_group_4fold
- weighted_top3_ensemble: Spearman 0.7147, RMSE 2.0205
- extratrees: Spearman 0.7138, RMSE 2.0263
- randomforest: Spearman 0.7111, RMSE 2.0197
- xgboost: Spearman 0.6740, RMSE 2.1296
- lightgbm: Spearman 0.6611, RMSE 2.1445

## legacy_rich_all_drugs_zero_smiles

- Input shape: [20353, 6574]
- Feature count: 5090
- Drugs: 295
- Cell lines: 76
- Invalid SMILES pairs kept: 7305

### drug_group_4fold
- weighted_top3_ensemble: Spearman 0.6396, RMSE 1.9928
- extratrees: Spearman 0.6256, RMSE 2.0071
- lightgbm: Spearman 0.6210, RMSE 2.0714
- xgboost: Spearman 0.6159, RMSE 2.0561
- randomforest: Spearman 0.6153, RMSE 2.0223

### random_4fold
- lightgbm: Spearman 0.8188, RMSE 1.3529
- weighted_top3_ensemble: Spearman 0.8114, RMSE 1.3948
- xgboost: Spearman 0.7920, RMSE 1.4560
- flat_mlp: Spearman 0.7864, RMSE 1.4777
- residual_mlp: Spearman 0.7852, RMSE 1.4948

### scaffold_group_4fold
- weighted_top3_ensemble: Spearman 0.6181, RMSE 2.0670
- randomforest: Spearman 0.6095, RMSE 2.0921
- extratrees: Spearman 0.6054, RMSE 2.0767
- lightgbm: Spearman 0.5878, RMSE 2.1327
- xgboost: Spearman 0.5664, RMSE 2.1773

## colonstyle_compact_baseline

- Input shape: [20353, 6574]
- Feature count: 4735
- Drugs: 295
- Cell lines: 76
- Invalid SMILES pairs kept: 7305

### drug_group_4fold
- weighted_top3_ensemble: Spearman 0.5661, RMSE 2.1594
- xgboost: Spearman 0.5488, RMSE 2.1706
- lightgbm: Spearman 0.5485, RMSE 2.1747
- randomforest: Spearman 0.5219, RMSE 2.2796
- extratrees: Spearman 0.5099, RMSE 2.2920

### random_4fold
- weighted_top3_ensemble: Spearman 0.7883, RMSE 1.4722
- lightgbm: Spearman 0.7862, RMSE 1.4823
- flat_mlp: Spearman 0.7756, RMSE 1.5064
- residual_mlp: Spearman 0.7712, RMSE 1.5373
- xgboost: Spearman 0.7649, RMSE 1.5538

### scaffold_group_4fold
- weighted_top3_ensemble: Spearman 0.5299, RMSE 2.2085
- lightgbm: Spearman 0.5190, RMSE 2.2073
- randomforest: Spearman 0.5060, RMSE 2.3244
- xgboost: Spearman 0.4887, RMSE 2.2684
- extratrees: Spearman 0.4752, RMSE 2.3630

## Variant Policy Comparison

| Split | Best valid-only | Best all-drugs-zero | Delta valid-zero |
|---|---:|---:|---:|
| drug_group_4fold | 0.7430 | 0.6396 | +0.1034 |
| random_4fold | 0.8978 | 0.8188 | +0.0790 |
| scaffold_group_4fold | 0.7147 | 0.6181 | +0.0966 |

## Recommended Input Policy

- Recommended variant: `legacy_rich_valid_smiles_only`
- Basis: highest mean Spearman across drug_group_4fold and scaffold_group_4fold

| Variant | Mean drug/scaffold Spearman |
|---|---:|
| legacy_rich_valid_smiles_only | 0.7288 |
| legacy_rich_all_drugs_zero_smiles | 0.6288 |
| colonstyle_compact_baseline | 0.5480 |

## Extended Metrics For Split Winners

| Split | Variant | Model | Spearman | Pearson | Kendall | RMSE | MAE | R2 | AUROC | AUPRC |
|---|---|---|---:|---:|---:|---:|---:|---:|---:|---:|
| drug_group_4fold | legacy_rich_valid_smiles_only | extratrees | 0.7430 | 0.7570 | 0.5552 | 1.9388 | 1.4612 | 0.5682 | 0.8795 | 0.8286 |
| random_4fold | legacy_rich_valid_smiles_only | weighted_top3_ensemble | 0.8978 | 0.9216 | 0.7303 | 1.1544 | 0.8744 | 0.8469 | 0.9514 | 0.9214 |
| scaffold_group_4fold | legacy_rich_valid_smiles_only | weighted_top3_ensemble | 0.7147 | 0.7373 | 0.5283 | 2.0205 | 1.5150 | 0.5310 | 0.8669 | 0.8106 |

Binary metrics use `label_binary=1` as sensitive and score `-predicted LN_IC50`.
