# BRCA Hybrid Pipeline Report

- Cancer: Breast invasive carcinoma (유방암)
- Run ID: `260423_BRCA_V2_raw_lincs_cmap2020_all`
- Raw S3: `s3://say2-4team/oringinal_raw`
- Upload S3: `s3://say2-4team/20260409_eseo/260423_BRCA_V2/raw_lincs_cmap2020_all`
- Fold policy: 4-fold random/drug/scaffold

## Input QC

- Labels: 13106 rows, 51 cell lines, 295 drugs
- Valid SMILES drugs: 186 / 295
- DepMap mapped cell lines: 51 / 51
- CRISPR features: 3500
- LINCS mapped drugs: 176

## LINCS Policy

- Strategy: `all_cell_main`
- Main-run eligible: no
- Disease-overlap LINCS cells: BT20, BT474, HS578T, MCF7, MDAMB231, MDAMB468, T47D
- Representative cells declared: MCF7
- Matched drugs: 176
- Matched signatures: 81207
- Reason: This run keeps all-cell LINCS although 7 overlapping disease LINCS cell lines were available (BT20, BT474, HS578T, MCF7, MDAMB231, MDAMB468, T47D). In that situation, a disease-restricted LINCS main run should usually be preferred.

## legacy_rich_valid_smiles_only

- Input shape: [8409, 6551]
- Feature count: 5064
- Drugs: 186
- Cell lines: 51
- Invalid SMILES pairs kept: 0

### drug_group_4fold
- weighted_top3_ensemble: Spearman 0.6893, RMSE 2.0752
- extratrees: Spearman 0.6687, RMSE 2.0993
- randomforest: Spearman 0.6529, RMSE 2.1455
- flat_mlp: Spearman 0.6444, RMSE 2.1213
- xgboost: Spearman 0.6412, RMSE 2.2129

### random_4fold
- lightgbm: Spearman 0.8707, RMSE 1.2227
- weighted_top3_ensemble: Spearman 0.8703, RMSE 1.2339
- xgboost: Spearman 0.8605, RMSE 1.2757
- flat_mlp: Spearman 0.8551, RMSE 1.3003
- randomforest: Spearman 0.8291, RMSE 1.4094

### scaffold_group_4fold
- weighted_top3_ensemble: Spearman 0.7036, RMSE 1.9691
- extratrees: Spearman 0.6749, RMSE 2.0038
- flat_mlp: Spearman 0.6721, RMSE 1.9888
- randomforest: Spearman 0.6656, RMSE 2.0717
- residual_mlp: Spearman 0.6421, RMSE 2.0491

## legacy_rich_all_drugs_zero_smiles

- Input shape: [13106, 6551]
- Feature count: 5067
- Drugs: 295
- Cell lines: 51
- Invalid SMILES pairs kept: 4697

### drug_group_4fold
- weighted_top3_ensemble: Spearman 0.5809, RMSE 2.0399
- extratrees: Spearman 0.5773, RMSE 2.0477
- randomforest: Spearman 0.5719, RMSE 2.0591
- lightgbm: Spearman 0.5427, RMSE 2.1421
- xgboost: Spearman 0.5405, RMSE 2.1219

### random_4fold
- lightgbm: Spearman 0.7893, RMSE 1.4080
- weighted_top3_ensemble: Spearman 0.7838, RMSE 1.4462
- xgboost: Spearman 0.7679, RMSE 1.4890
- flat_mlp: Spearman 0.7559, RMSE 1.5390
- residual_mlp: Spearman 0.7349, RMSE 1.6075

### scaffold_group_4fold
- weighted_top3_ensemble: Spearman 0.5591, RMSE 2.0538
- xgboost: Spearman 0.5422, RMSE 2.0935
- randomforest: Spearman 0.5416, RMSE 2.0879
- lightgbm: Spearman 0.5402, RMSE 2.1012
- extratrees: Spearman 0.5309, RMSE 2.1018

## colonstyle_compact_baseline

- Input shape: [13106, 6551]
- Feature count: 4712
- Drugs: 295
- Cell lines: 51
- Invalid SMILES pairs kept: 4697

### drug_group_4fold
- weighted_top3_ensemble: Spearman 0.5119, RMSE 2.1513
- xgboost: Spearman 0.4968, RMSE 2.1989
- randomforest: Spearman 0.4871, RMSE 2.2014
- lightgbm: Spearman 0.4813, RMSE 2.2191
- extratrees: Spearman 0.4715, RMSE 2.2404

### random_4fold
- weighted_top3_ensemble: Spearman 0.7546, RMSE 1.5338
- lightgbm: Spearman 0.7533, RMSE 1.5324
- xgboost: Spearman 0.7401, RMSE 1.5784
- flat_mlp: Spearman 0.7367, RMSE 1.5667
- residual_mlp: Spearman 0.7287, RMSE 1.6000

### scaffold_group_4fold
- weighted_top3_ensemble: Spearman 0.4950, RMSE 2.1594
- xgboost: Spearman 0.4947, RMSE 2.1810
- lightgbm: Spearman 0.4637, RMSE 2.1919
- randomforest: Spearman 0.4558, RMSE 2.2675
- extratrees: Spearman 0.4369, RMSE 2.2960

## Variant Policy Comparison

| Split | Best valid-only | Best all-drugs-zero | Delta valid-zero |
|---|---:|---:|---:|
| drug_group_4fold | 0.6893 | 0.5809 | +0.1084 |
| random_4fold | 0.8707 | 0.7893 | +0.0814 |
| scaffold_group_4fold | 0.7036 | 0.5591 | +0.1445 |

## Recommended Input Policy

- Recommended variant: `legacy_rich_valid_smiles_only`
- Basis: highest mean Spearman across drug_group_4fold and scaffold_group_4fold

| Variant | Mean drug/scaffold Spearman |
|---|---:|
| legacy_rich_valid_smiles_only | 0.6965 |
| legacy_rich_all_drugs_zero_smiles | 0.5700 |
| colonstyle_compact_baseline | 0.5034 |

## Extended Metrics For Split Winners

| Split | Variant | Model | Spearman | Pearson | Kendall | RMSE | MAE | R2 | AUROC | AUPRC |
|---|---|---|---:|---:|---:|---:|---:|---:|---:|---:|
| drug_group_4fold | legacy_rich_valid_smiles_only | weighted_top3_ensemble | 0.6893 | 0.6944 | 0.5034 | 2.0752 | 1.5467 | 0.4732 | 0.8540 | 0.7817 |
| random_4fold | legacy_rich_valid_smiles_only | lightgbm | 0.8707 | 0.9044 | 0.6971 | 1.2227 | 0.9094 | 0.8171 | 0.9406 | 0.9018 |
| scaffold_group_4fold | legacy_rich_valid_smiles_only | weighted_top3_ensemble | 0.7036 | 0.7413 | 0.5214 | 1.9691 | 1.4941 | 0.5257 | 0.8547 | 0.7916 |

Binary metrics use `label_binary=1` as sensitive and score `-predicted LN_IC50`.
