# Liver LINCS Policy Comparison

- Date: `2026-04-23`
- Main comparison runs:
  - `260423_LIVER_V2_raw_lincs_cmap2020_all`
  - `260423_LIVER_V2_raw_lincs_cmap2020_overlap`

## Structural policy tags

- all-cell run: `all_cell_main`
- overlap run: `cancer_specific_main`

## Coverage

### All-cell

- LINCS mapped drugs: `176`
- matched signatures: `81207`
- included LINCS cell lines: `221`

### Direct overlap (`HUH7`, `JHH7`)

- LINCS mapped drugs: `61`
- matched signatures: `184`
- included LINCS cell lines: `2`

## Best split winners

| Split | All-cell best | Overlap best | Delta |
|---|---:|---:|---:|
| random_4fold | 0.8752 | 0.8680 | +0.0072 |
| drug_group_4fold | 0.6910 | 0.5590 | +0.1320 |
| scaffold_group_4fold | 0.6712 | 0.5435 | +0.1278 |

## Variant recommendation

Both runs recommended:

- `legacy_rich_valid_smiles_only`

## Conclusion

For liver, the final practical direction is:

- `all-cell LINCS`
- `legacy_rich_valid_smiles_only`

Even though a direct LIHC overlap run exists (`HUH7`, `JHH7`), it is materially weaker than the all-cell run on drug/scaffold holdout and should remain a supporting analysis rather than the main setting.
