# Calculate allele prevalence

Calculate allele prevalence

## Usage

``` r
calculate_prevalence(allele_table)
```

## Arguments

- allele_table:

  Allele table with `specimen_name`, `target_name`, and `variant`.

## Value

A tibble with `target_name`, `variant`, `prev`, `sample_count`, and
`sample_total`.
