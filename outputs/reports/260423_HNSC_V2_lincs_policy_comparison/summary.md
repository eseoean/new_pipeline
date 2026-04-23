# HNSC LINCS Policy Comparison

- Date: `2026-04-23`
- Main run: `260423_HNSC_V2_raw_lincs_cmap2020_all`
- Supporting run: `260423_HNSC_V2_raw_lincs_cmap2020_bicr6`

## Policy Outcome

- all-cell run: `all_cell_main`
- BICR6 explicit run: `cancer_specific_unavailable`

## Coverage

### All-cell

- LINCS mapped drugs: `176`
- matched signatures: `81207`
- included LINCS cell lines: `221`
- disease-overlap LINCS cells: `0`

### BICR6 explicit

- LINCS mapped drugs: `38`
- matched signatures: `114`
- included LINCS cell lines: `1`

## Best split winners

| Split | All-cell best | BICR6 best | Delta |
|---|---:|---:|---:|
| random_4fold | 0.9091 | 0.9056 | +0.0035 |
| drug_group_4fold | 0.7317 | 0.5720 | +0.1597 |
| scaffold_group_4fold | 0.7063 | 0.5765 | +0.1297 |

## Variant recommendation

Both runs recommended:

- `legacy_rich_valid_smiles_only`

## Conclusion

For HNSC, the current adaptive LINCS rule selects:

- `all_cell_main`

`BICR6` is biologically suggestive, but it is not a direct training-cell overlap and its coverage is too small to replace the stronger all-cell LINCS block.
