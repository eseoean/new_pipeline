# BRCA Hybrid Pipeline Report

- Cancer: Breast invasive carcinoma (유방암)
- Run ID: `260423_BRCA_V2_raw_lincs_cmap2020_overlap`
- Raw S3: `s3://say2-4team/oringinal_raw`
- Upload S3: `s3://say2-4team/20260409_eseo/260423_BRCA_V2/raw_lincs_cmap2020_overlap`
- Fold policy: 4-fold random/drug/scaffold

## Input QC

- Labels: 13106 rows, 51 cell lines, 295 drugs
- Valid SMILES drugs: 186 / 295
- DepMap mapped cell lines: 51 / 51
- CRISPR features: 3500
- LINCS mapped drugs: 152

## LINCS Policy

- Strategy: `cancer_specific_main`
- Main-run eligible: yes
- Disease-overlap LINCS cells: BT20, BT474, HS578T, MCF7, MDAMB231, MDAMB468, T47D
- Representative cells declared: MCF7
- Matched drugs: 152
- Matched signatures: 11340
- Reason: 7 overlapping disease LINCS cell lines were available (BT20, BT474, HS578T, MCF7, MDAMB231, MDAMB468, T47D), so cancer-specific LINCS can be used as the main setting.

## legacy_rich_valid_smiles_only

- Input shape: [8409, 6551]
- Feature count: 5064
- Drugs: 186
- Cell lines: 51
- Invalid SMILES pairs kept: 0

### drug_group_4fold
- weighted_top3_ensemble: Spearman 0.6556, RMSE 2.1478
- randomforest: Spearman 0.6388, RMSE 2.1654
- extratrees: Spearman 0.6377, RMSE 2.1497
- xgboost: Spearman 0.6226, RMSE 2.2848
- lightgbm: Spearman 0.5876, RMSE 2.3900

### random_4fold
- weighted_top3_ensemble: Spearman 0.8692, RMSE 1.2407
- lightgbm: Spearman 0.8680, RMSE 1.2344
- xgboost: Spearman 0.8590, RMSE 1.2794
- flat_mlp: Spearman 0.8541, RMSE 1.3154
- randomforest: Spearman 0.8268, RMSE 1.4174

### scaffold_group_4fold
- weighted_top3_ensemble: Spearman 0.6679, RMSE 2.0379
- extratrees: Spearman 0.6530, RMSE 2.0451
- xgboost: Spearman 0.6270, RMSE 2.1751
- flat_mlp: Spearman 0.6266, RMSE 2.0763
- randomforest: Spearman 0.6223, RMSE 2.1175

## legacy_rich_all_drugs_zero_smiles

- Input shape: [13106, 6551]
- Feature count: 5067
- Drugs: 295
- Cell lines: 51
- Invalid SMILES pairs kept: 4697

### drug_group_4fold
- weighted_top3_ensemble: Spearman 0.5541, RMSE 2.0723
- extratrees: Spearman 0.5356, RMSE 2.0824
- xgboost: Spearman 0.5351, RMSE 2.1676
- randomforest: Spearman 0.5335, RMSE 2.1086
- lightgbm: Spearman 0.5229, RMSE 2.1918

### random_4fold
- lightgbm: Spearman 0.7828, RMSE 1.4208
- weighted_top3_ensemble: Spearman 0.7777, RMSE 1.4577
- xgboost: Spearman 0.7593, RMSE 1.5063
- flat_mlp: Spearman 0.7486, RMSE 1.5541
- residual_mlp: Spearman 0.7354, RMSE 1.6066

### scaffold_group_4fold
- weighted_top3_ensemble: Spearman 0.5395, RMSE 2.1025
- xgboost: Spearman 0.5339, RMSE 2.1554
- lightgbm: Spearman 0.5195, RMSE 2.1910
- extratrees: Spearman 0.5086, RMSE 2.1292
- randomforest: Spearman 0.5017, RMSE 2.1555

## colonstyle_compact_baseline

- Input shape: [13106, 6551]
- Feature count: 4712
- Drugs: 295
- Cell lines: 51
- Invalid SMILES pairs kept: 4697

### drug_group_4fold
- weighted_top3_ensemble: Spearman 0.5139, RMSE 2.1538
- lightgbm: Spearman 0.5029, RMSE 2.2092
- xgboost: Spearman 0.4915, RMSE 2.1985
- randomforest: Spearman 0.4841, RMSE 2.2110
- extratrees: Spearman 0.4775, RMSE 2.2330

### random_4fold
- weighted_top3_ensemble: Spearman 0.7517, RMSE 1.5382
- lightgbm: Spearman 0.7488, RMSE 1.5406
- xgboost: Spearman 0.7356, RMSE 1.5864
- flat_mlp: Spearman 0.7352, RMSE 1.5669
- residual_mlp: Spearman 0.7312, RMSE 1.5950

### scaffold_group_4fold
- weighted_top3_ensemble: Spearman 0.4831, RMSE 2.2002
- lightgbm: Spearman 0.4717, RMSE 2.2371
- xgboost: Spearman 0.4711, RMSE 2.2549
- randomforest: Spearman 0.4446, RMSE 2.2688
- extratrees: Spearman 0.4288, RMSE 2.3075

## Variant Policy Comparison

| Split | Best valid-only | Best all-drugs-zero | Delta valid-zero |
|---|---:|---:|---:|
| drug_group_4fold | 0.6556 | 0.5541 | +0.1014 |
| random_4fold | 0.8692 | 0.7828 | +0.0864 |
| scaffold_group_4fold | 0.6679 | 0.5395 | +0.1284 |

## Recommended Input Policy

- Recommended variant: `legacy_rich_valid_smiles_only`
- Basis: highest mean Spearman across drug_group_4fold and scaffold_group_4fold

| Variant | Mean drug/scaffold Spearman |
|---|---:|
| legacy_rich_valid_smiles_only | 0.6617 |
| legacy_rich_all_drugs_zero_smiles | 0.5468 |
| colonstyle_compact_baseline | 0.4985 |

## Extended Metrics For Split Winners

| Split | Variant | Model | Spearman | Pearson | Kendall | RMSE | MAE | R2 | AUROC | AUPRC |
|---|---|---|---:|---:|---:|---:|---:|---:|---:|---:|
| drug_group_4fold | legacy_rich_valid_smiles_only | weighted_top3_ensemble | 0.6556 | 0.6641 | 0.4778 | 2.1478 | 1.6207 | 0.4357 | 0.8299 | 0.7414 |
| random_4fold | legacy_rich_valid_smiles_only | weighted_top3_ensemble | 0.8692 | 0.9029 | 0.6949 | 1.2407 | 0.9237 | 0.8117 | 0.9400 | 0.9009 |
| scaffold_group_4fold | legacy_rich_valid_smiles_only | weighted_top3_ensemble | 0.6679 | 0.7150 | 0.4915 | 2.0379 | 1.5644 | 0.4920 | 0.8331 | 0.7603 |

Binary metrics use `label_binary=1` as sensitive and score `-predicted LN_IC50`.
