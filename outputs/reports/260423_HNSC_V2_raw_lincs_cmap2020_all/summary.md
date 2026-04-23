# HNSC Hybrid Pipeline Report

- Cancer: Head and neck squamous cell carcinoma (두경부 편평세포암)
- Run ID: `260423_HNSC_V2_raw_lincs_cmap2020_all`
- Raw S3: `s3://say2-4team/HNSC_raw`
- Upload S3: `s3://say2-4team/20260409_eseo/260423_HNSC_V2/raw_lincs_cmap2020_all`
- Fold policy: 4-fold random/drug/scaffold

## Input QC

- Labels: 9358 rows, 39 cell lines, 295 drugs
- Valid SMILES drugs: 186 / 295
- DepMap mapped cell lines: 39 / 39
- CRISPR features: 3500
- LINCS mapped drugs: 176

## LINCS Policy

- Strategy: `all_cell_main`
- Main-run eligible: yes
- Disease-overlap LINCS cells: none
- Representative cells declared: none
- Matched drugs: 176
- Matched signatures: 81207
- Reason: No directly overlapping disease LINCS cell line was available after drug matching, so all-cell drug-level LINCS is used as the fallback main setting.

## legacy_rich_valid_smiles_only

- Input shape: [6259, 6539]
- Feature count: 5052
- Drugs: 186
- Cell lines: 39
- Invalid SMILES pairs kept: 0

### drug_group_4fold
- weighted_top3_ensemble: Spearman 0.7317, RMSE 1.9368
- extratrees: Spearman 0.7288, RMSE 1.9410
- randomforest: Spearman 0.7221, RMSE 1.9433
- xgboost: Spearman 0.6996, RMSE 2.0513
- flat_mlp: Spearman 0.6936, RMSE 2.0012

### random_4fold
- weighted_top3_ensemble: Spearman 0.9091, RMSE 1.0305
- lightgbm: Spearman 0.9079, RMSE 1.0283
- xgboost: Spearman 0.9041, RMSE 1.0562
- flat_mlp: Spearman 0.8977, RMSE 1.1003
- randomforest: Spearman 0.8910, RMSE 1.1451

### scaffold_group_4fold
- weighted_top3_ensemble: Spearman 0.7063, RMSE 2.1093
- extratrees: Spearman 0.6944, RMSE 2.1241
- randomforest: Spearman 0.6906, RMSE 2.1329
- xgboost: Spearman 0.6822, RMSE 2.1882
- lightgbm: Spearman 0.6734, RMSE 2.1951

## legacy_rich_all_drugs_zero_smiles

- Input shape: [9358, 6539]
- Feature count: 5055
- Drugs: 295
- Cell lines: 39
- Invalid SMILES pairs kept: 3099

### drug_group_4fold
- weighted_top3_ensemble: Spearman 0.6268, RMSE 1.9425
- extratrees: Spearman 0.6230, RMSE 1.9451
- randomforest: Spearman 0.6178, RMSE 1.9631
- xgboost: Spearman 0.5929, RMSE 2.0538
- flat_mlp: Spearman 0.5736, RMSE 1.9722

### random_4fold
- lightgbm: Spearman 0.8264, RMSE 1.2635
- weighted_top3_ensemble: Spearman 0.8213, RMSE 1.3008
- xgboost: Spearman 0.8125, RMSE 1.3219
- flat_mlp: Spearman 0.7911, RMSE 1.4056
- residual_mlp: Spearman 0.7823, RMSE 1.4515

### scaffold_group_4fold
- weighted_top3_ensemble: Spearman 0.6237, RMSE 1.9376
- extratrees: Spearman 0.6116, RMSE 1.9740
- randomforest: Spearman 0.5943, RMSE 2.0314
- flat_mlp: Spearman 0.5859, RMSE 1.9691
- residual_mlp: Spearman 0.5839, RMSE 2.0113

## colonstyle_compact_baseline

- Input shape: [9358, 6539]
- Feature count: 4700
- Drugs: 295
- Cell lines: 39
- Invalid SMILES pairs kept: 3099

### drug_group_4fold
- weighted_top3_ensemble: Spearman 0.5500, RMSE 2.0605
- randomforest: Spearman 0.5310, RMSE 2.1703
- xgboost: Spearman 0.5253, RMSE 2.1040
- lightgbm: Spearman 0.5185, RMSE 2.1189
- extratrees: Spearman 0.4927, RMSE 2.2206

### random_4fold
- lightgbm: Spearman 0.7954, RMSE 1.3804
- weighted_top3_ensemble: Spearman 0.7896, RMSE 1.4003
- xgboost: Spearman 0.7800, RMSE 1.4322
- residual_mlp: Spearman 0.7602, RMSE 1.4862
- flat_mlp: Spearman 0.7593, RMSE 1.4640

### scaffold_group_4fold
- weighted_top3_ensemble: Spearman 0.5084, RMSE 2.2322
- randomforest: Spearman 0.4887, RMSE 2.2965
- extratrees: Spearman 0.4733, RMSE 2.3101
- xgboost: Spearman 0.4633, RMSE 2.3058
- lightgbm: Spearman 0.4379, RMSE 2.3786

## Variant Policy Comparison

| Split | Best valid-only | Best all-drugs-zero | Delta valid-zero |
|---|---:|---:|---:|
| drug_group_4fold | 0.7317 | 0.6268 | +0.1049 |
| random_4fold | 0.9091 | 0.8264 | +0.0827 |
| scaffold_group_4fold | 0.7063 | 0.6237 | +0.0826 |

## Recommended Input Policy

- Recommended variant: `legacy_rich_valid_smiles_only`
- Basis: highest mean Spearman across drug_group_4fold and scaffold_group_4fold

| Variant | Mean drug/scaffold Spearman |
|---|---:|
| legacy_rich_valid_smiles_only | 0.7190 |
| legacy_rich_all_drugs_zero_smiles | 0.6253 |
| colonstyle_compact_baseline | 0.5292 |

## Extended Metrics For Split Winners

| Split | Variant | Model | Spearman | Pearson | Kendall | RMSE | MAE | R2 | AUROC | AUPRC |
|---|---|---|---:|---:|---:|---:|---:|---:|---:|---:|
| drug_group_4fold | legacy_rich_valid_smiles_only | weighted_top3_ensemble | 0.7317 | 0.7444 | 0.5416 | 1.9368 | 1.4371 | 0.5505 | 0.8746 | 0.8128 |
| random_4fold | legacy_rich_valid_smiles_only | weighted_top3_ensemble | 0.9091 | 0.9349 | 0.7486 | 1.0305 | 0.7647 | 0.8728 | 0.9597 | 0.9340 |
| scaffold_group_4fold | legacy_rich_valid_smiles_only | weighted_top3_ensemble | 0.7063 | 0.6941 | 0.5143 | 2.1093 | 1.5163 | 0.4669 | 0.8610 | 0.7877 |

Binary metrics use `label_binary=1` as sensitive and score `-predicted LN_IC50`.
