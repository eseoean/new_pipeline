# Lung LINCS Policy Comparison

- Date: `2026-04-23`
- Main comparison runs:
  - `260423_LUNG_V2_raw_lincs_cmap2020_all`
  - `260423_LUNG_V2_raw_lincs_cmap2020_overlap`
- Disease scope: `LUAD + LUSC only`
- Excluded from label scope: `SCLC`

## Structural policy tags

- all-cell run: `all_cell_main`
- overlap run: `cancer_specific_main`

## Coverage

### All-cell

- LINCS mapped drugs: `176`
- matched signatures: `81207`
- included LINCS cell lines: `221`

### Direct NSCLC overlap (`10` cell lines)

- LINCS mapped drugs: `152`
- matched signatures: `5649`
- included LINCS cell lines: `10`

## Best split winners

| Split | All-cell best | Overlap best | Delta |
|---|---:|---:|---:|
| random_4fold | 0.8978 | 0.8976 | +0.0002 |
| drug_group_4fold | 0.7430 | 0.7026 | +0.0404 |
| scaffold_group_4fold | 0.7147 | 0.7161 | -0.0015 |

## Variant recommendation

Both runs recommended:

- `legacy_rich_valid_smiles_only`

## Conclusion

For lung NSCLC, the final practical direction is:

- `all-cell LINCS`
- `legacy_rich_valid_smiles_only`

This is a good example of why the adaptive policy should be benchmark-driven at the end:

- structurally, a strong disease-overlap run exists
- empirically, all-cell still gives the better overall generalization profile because the drug-holdout gain is larger than the tiny scaffold-holdout loss
