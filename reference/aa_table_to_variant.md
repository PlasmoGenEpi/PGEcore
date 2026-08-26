# Convert long-form amino acid calls to variant strings

Convert long-form amino acid calls to variant strings

## Usage

``` r
aa_table_to_variant(aa_table)
```

## Arguments

- aa_table:

  Tibble with `specimen_name`, `gene`, `pos`, `reads`, `aa`, `n_aa`.

## Value

Character vector of variant strings, one per specimen.
