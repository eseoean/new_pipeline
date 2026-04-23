# LIHC Hybrid Pipeline Report

- Cancer: Liver hepatocellular carcinoma (간암)
- Run ID: `260423_LIVER_V2_raw_lincs_cmap2020_overlap`
- Raw S3: `s3://say2-4team/Liver_raw`
- Upload S3: `s3://say2-4team/20260409_eseo/260423_LIVER_V2/raw_lincs_cmap2020_overlap`
- Fold policy: 4-fold random/drug/scaffold

## Input QC

- Labels: 4164 rows, 15 cell lines, 295 drugs
- Valid SMILES drugs: 186 / 295
- DepMap mapped cell lines: 15 / 15
- CRISPR features: 3500
- LINCS mapped drugs: 61

## LINCS Policy

- Strategy: `cancer_specific_main`
- Main-run eligible: yes
- Disease-overlap LINCS cells: HUH7, JHH7
- Representative cells declared: none
- Matched drugs: 61
- Matched signatures: 184
- Reason: 2 overlapping disease LINCS cell lines were available (HUH7, JHH7), so cancer-specific LINCS can be used as the main setting.

## legacy_rich_valid_smiles_only

- Input shape: [2590, 6511]
- Feature count: 5024
- Drugs: 186
- Cell lines: 15
- Invalid SMILES pairs kept: 0

### drug_group_4fold
- weighted_top3_ensemble: Spearman 0.5590, RMSE 2.2049
- randomforest: Spearman 0.5391, RMSE 2.2214
- extratrees: Spearman 0.5298, RMSE 2.2383
- lightgbm: Spearman 0.5184, RMSE 2.3105
- xgboost: Spearman 0.4933, RMSE 2.3445

### random_4fold
- lightgbm: Spearman 0.8680, RMSE 1.1956
- weighted_top3_ensemble: Spearman 0.8644, RMSE 1.2129
- xgboost: Spearman 0.8619, RMSE 1.2175
- randomforest: Spearman 0.8343, RMSE 1.3179
- extratrees: Spearman 0.8332, RMSE 1.3405

### scaffold_group_4fold
- weighted_top3_ensemble: Spearman 0.5435, RMSE 2.3202
- randomforest: Spearman 0.5127, RMSE 2.3597
- residual_mlp: Spearman 0.4680, RMSE 2.3969
- xgboost: Spearman 0.4674, RMSE 2.4314
- extratrees: Spearman 0.4636, RMSE 2.3906

## legacy_rich_all_drugs_zero_smiles

- Input shape: [4164, 6511]
- Feature count: 5027
- Drugs: 295
- Cell lines: 15
- Invalid SMILES pairs kept: 1574

### drug_group_4fold
- weighted_top3_ensemble: Spearman 0.4349, RMSE 2.2648
- xgboost: Spearman 0.4182, RMSE 2.3041
- lightgbm: Spearman 0.4126, RMSE 2.3327
- extratrees: Spearman 0.4000, RMSE 2.2922
- randomforest: Spearman 0.3958, RMSE 2.2867

### random_4fold
- lightgbm: Spearman 0.7594, RMSE 1.4352
- xgboost: Spearman 0.7509, RMSE 1.4783
- weighted_top3_ensemble: Spearman 0.7460, RMSE 1.4803
- randomforest: Spearman 0.6770, RMSE 1.6334
- extratrees: Spearman 0.6643, RMSE 1.6643

### scaffold_group_4fold
- weighted_top3_ensemble: Spearman 0.4500, RMSE 2.2360
- randomforest: Spearman 0.4445, RMSE 2.2474
- extratrees: Spearman 0.4221, RMSE 2.2742
- xgboost: Spearman 0.4133, RMSE 2.3025
- lightgbm: Spearman 0.4002, RMSE 2.3618

## colonstyle_compact_baseline

- Input shape: [4164, 6511]
- Feature count: 4672
- Drugs: 295
- Cell lines: 15
- Invalid SMILES pairs kept: 1574

### drug_group_4fold
- extratrees: Spearman 0.3846, RMSE 2.3298
- weighted_top3_ensemble: Spearman 0.3798, RMSE 2.3213
- randomforest: Spearman 0.3635, RMSE 2.3380
- xgboost: Spearman 0.3470, RMSE 2.4095
- residual_mlp: Spearman 0.3295, RMSE 2.3608

### random_4fold
- lightgbm: Spearman 0.7268, RMSE 1.5494
- weighted_top3_ensemble: Spearman 0.7239, RMSE 1.5620
- xgboost: Spearman 0.7183, RMSE 1.5728
- flat_mlp: Spearman 0.6602, RMSE 1.6802
- randomforest: Spearman 0.6471, RMSE 1.7160

### scaffold_group_4fold
- weighted_top3_ensemble: Spearman 0.3766, RMSE 2.3277
- xgboost: Spearman 0.3748, RMSE 2.3454
- randomforest: Spearman 0.3592, RMSE 2.3541
- extratrees: Spearman 0.3317, RMSE 2.3854
- lightgbm: Spearman 0.3149, RMSE 2.4713

## Variant Policy Comparison

| Split | Best valid-only | Best all-drugs-zero | Delta valid-zero |
|---|---:|---:|---:|
| drug_group_4fold | 0.5590 | 0.4349 | +0.1241 |
| random_4fold | 0.8680 | 0.7594 | +0.1086 |
| scaffold_group_4fold | 0.5435 | 0.4500 | +0.0934 |

## Recommended Input Policy

- Recommended variant: `legacy_rich_valid_smiles_only`
- Basis: highest mean Spearman across drug_group_4fold and scaffold_group_4fold

| Variant | Mean drug/scaffold Spearman |
|---|---:|
| legacy_rich_valid_smiles_only | 0.5512 |
| legacy_rich_all_drugs_zero_smiles | 0.4425 |
| colonstyle_compact_baseline | 0.3806 |

## Extended Metrics For Split Winners

| Split | Variant | Model | Spearman | Pearson | Kendall | RMSE | MAE | R2 | AUROC | AUPRC |
|---|---|---|---:|---:|---:|---:|---:|---:|---:|---:|
| drug_group_4fold | legacy_rich_valid_smiles_only | weighted_top3_ensemble | 0.5590 | 0.5947 | 0.3965 | 2.2049 | 1.6739 | 0.3531 | 0.7816 | 0.6895 |
| random_4fold | legacy_rich_valid_smiles_only | lightgbm | 0.8680 | 0.9000 | 0.6909 | 1.1956 | 0.9118 | 0.8098 | 0.9428 | 0.8991 |
| scaffold_group_4fold | legacy_rich_valid_smiles_only | weighted_top3_ensemble | 0.5435 | 0.5477 | 0.3846 | 2.3202 | 1.7242 | 0.2837 | 0.7806 | 0.6720 |

Binary metrics use `label_binary=1` as sensitive and score `-predicted LN_IC50`.
