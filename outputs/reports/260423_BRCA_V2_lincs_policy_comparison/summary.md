# BRCA LINCS Policy Comparison

- Date: `2026-04-23`
- Main comparison runs:
  - `260423_BRCA_V2_raw_lincs_cmap2020_all`
  - `260423_BRCA_V2_raw_lincs_cmap2020_overlap`
  - `260423_BRCA_V2_raw_lincs_cmap2020_mcf7`
- Disease scope: `BRCA`
- Representative benchmark: `MCF7`

## Structural policy tags

- all-cell run: `all_cell_main`
- overlap run: `cancer_specific_main`
- MCF7 explicit run: representative-line benchmark

## Coverage

### All-cell

- LINCS mapped drugs: `176`
- matched signatures: `81207`
- included LINCS cell lines: `221`

### Direct BRCA overlap (`7` cell lines)

- LINCS mapped drugs: `152`
- matched signatures: `11340`
- included LINCS cell lines: `7`

### MCF7-only

- LINCS mapped drugs: `148`
- matched signatures: `6224`
- included LINCS cell lines: `1`

## Best split winners

| Split | All-cell best | Overlap best | MCF7-only best |
|---|---:|---:|---:|
| random_4fold | 0.8707 | 0.8692 | 0.8699 |
| drug_group_4fold | 0.6893 | 0.6556 | 0.6708 |
| scaffold_group_4fold | 0.7036 | 0.6679 | 0.6821 |

## Variant recommendation

All three runs recommended:

- `legacy_rich_valid_smiles_only`

## Conclusion

For BRCA, the final practical direction is:

- `all-cell LINCS`
- `legacy_rich_valid_smiles_only`

This is important because BRCA had the strongest reason so far to remain disease-specific or `MCF7`-only.
Even so, the benchmark still favored `all-cell` on the more important held-out drug and scaffold views.
