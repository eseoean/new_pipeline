# PAAD Hybrid Pipeline Report

- Cancer: Pancreatic adenocarcinoma (췌장암)
- Run ID: `260423_PAAD_V2_raw_lincs_cmap2020_yapc`
- Raw S3: `s3://say2-4team/PAAD_raw`
- Upload S3: `s3://say2-4team/20260409_eseo/260423_PAAD_V2/raw_lincs_cmap2020_yapc`
- Fold policy: 4-fold random/drug/scaffold

## Input QC

- Labels: 7513 rows, 29 cell lines, 295 drugs
- Valid SMILES drugs: 186 / 295
- DepMap mapped cell lines: 29 / 29
- CRISPR features: 3500
- LINCS mapped drugs: 116

## legacy_rich_valid_smiles_only

- Input shape: [4759, 6524]
- Feature count: 5036
- Drugs: 186
- Cell lines: 29
- Invalid SMILES pairs kept: 0

### drug_group_4fold
- weighted_top3_ensemble: Spearman 0.6887, RMSE 1.9415
- randomforest: Spearman 0.6782, RMSE 1.9788
- extratrees: Spearman 0.6666, RMSE 1.9926
- lightgbm: Spearman 0.6597, RMSE 2.0274
- xgboost: Spearman 0.6595, RMSE 2.0557

### random_4fold
- weighted_top3_ensemble: Spearman 0.9052, RMSE 1.0538
- lightgbm: Spearman 0.9034, RMSE 1.0581
- xgboost: Spearman 0.8994, RMSE 1.0888
- flat_mlp: Spearman 0.8929, RMSE 1.1129
- randomforest: Spearman 0.8924, RMSE 1.1456

### scaffold_group_4fold
- weighted_top3_ensemble: Spearman 0.6530, RMSE 2.1359
- randomforest: Spearman 0.6384, RMSE 2.1755
- extratrees: Spearman 0.6267, RMSE 2.1595
- xgboost: Spearman 0.6130, RMSE 2.2565
- flat_mlp: Spearman 0.6021, RMSE 2.1669

## legacy_rich_all_drugs_zero_smiles

- Input shape: [7513, 6524]
- Feature count: 5039
- Drugs: 295
- Cell lines: 29
- Invalid SMILES pairs kept: 2754

### drug_group_4fold
- weighted_top3_ensemble: Spearman 0.6100, RMSE 1.9757
- xgboost: Spearman 0.6043, RMSE 2.0025
- lightgbm: Spearman 0.5935, RMSE 2.0373
- randomforest: Spearman 0.5863, RMSE 1.9941
- extratrees: Spearman 0.5773, RMSE 2.0255

### random_4fold
- lightgbm: Spearman 0.8133, RMSE 1.2994
- weighted_top3_ensemble: Spearman 0.7987, RMSE 1.3545
- xgboost: Spearman 0.7935, RMSE 1.3719
- randomforest: Spearman 0.7522, RMSE 1.4903
- extratrees: Spearman 0.7420, RMSE 1.5289

### scaffold_group_4fold
- weighted_top3_ensemble: Spearman 0.5735, RMSE 2.0660
- randomforest: Spearman 0.5718, RMSE 2.0718
- extratrees: Spearman 0.5560, RMSE 2.1022
- xgboost: Spearman 0.5358, RMSE 2.1659
- lightgbm: Spearman 0.5153, RMSE 2.2500

## colonstyle_compact_baseline

- Input shape: [7513, 6524]
- Feature count: 4684
- Drugs: 295
- Cell lines: 29
- Invalid SMILES pairs kept: 2754

### drug_group_4fold
- weighted_top3_ensemble: Spearman 0.5870, RMSE 2.0546
- xgboost: Spearman 0.5764, RMSE 2.0426
- lightgbm: Spearman 0.5673, RMSE 2.0890
- randomforest: Spearman 0.5388, RMSE 2.1985
- extratrees: Spearman 0.5080, RMSE 2.2383

### random_4fold
- lightgbm: Spearman 0.7780, RMSE 1.4348
- weighted_top3_ensemble: Spearman 0.7676, RMSE 1.4644
- xgboost: Spearman 0.7671, RMSE 1.4644
- randomforest: Spearman 0.7315, RMSE 1.5705
- extratrees: Spearman 0.7234, RMSE 1.6026

### scaffold_group_4fold
- weighted_top3_ensemble: Spearman 0.5142, RMSE 2.1718
- xgboost: Spearman 0.5070, RMSE 2.1974
- randomforest: Spearman 0.4955, RMSE 2.2601
- lightgbm: Spearman 0.4855, RMSE 2.2506
- extratrees: Spearman 0.4666, RMSE 2.3205

## Variant Policy Comparison

| Split | Best valid-only | Best all-drugs-zero | Delta valid-zero |
|---|---:|---:|---:|
| drug_group_4fold | 0.6887 | 0.6100 | +0.0787 |
| random_4fold | 0.9052 | 0.8133 | +0.0919 |
| scaffold_group_4fold | 0.6530 | 0.5735 | +0.0794 |

## Recommended Input Policy

- Recommended variant: `legacy_rich_valid_smiles_only`
- Basis: highest mean Spearman across drug_group_4fold and scaffold_group_4fold

| Variant | Mean drug/scaffold Spearman |
|---|---:|
| legacy_rich_valid_smiles_only | 0.6708 |
| legacy_rich_all_drugs_zero_smiles | 0.5917 |
| colonstyle_compact_baseline | 0.5506 |

## Extended Metrics For Split Winners

| Split | Variant | Model | Spearman | Pearson | Kendall | RMSE | MAE | R2 | AUROC | AUPRC |
|---|---|---|---:|---:|---:|---:|---:|---:|---:|---:|
| drug_group_4fold | legacy_rich_valid_smiles_only | weighted_top3_ensemble | 0.6887 | 0.7472 | 0.5109 | 1.9415 | 1.5054 | 0.5529 | 0.8549 | 0.8050 |
| random_4fold | legacy_rich_valid_smiles_only | weighted_top3_ensemble | 0.9052 | 0.9330 | 0.7435 | 1.0538 | 0.8115 | 0.8683 | 0.9615 | 0.9359 |
| scaffold_group_4fold | legacy_rich_valid_smiles_only | weighted_top3_ensemble | 0.6530 | 0.6877 | 0.4762 | 2.1359 | 1.6205 | 0.4588 | 0.8299 | 0.7576 |

Binary metrics use `label_binary=1` as sensitive and score `-predicted LN_IC50`.
