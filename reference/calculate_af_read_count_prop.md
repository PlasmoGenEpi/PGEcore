# Calculate allele frequency from within-sample allele proportions

Calculate allele frequency from within-sample allele proportions

## Usage

``` r
calculate_af_read_count_prop(allele_table)
```

## Arguments

- allele_table:

  Allele table with `specimen_name`, `target_name`, `variant`, and
  `reads`.

## Value

A tibble with `target_name`, `variant`, and `freq`.
