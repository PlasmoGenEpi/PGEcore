# Calculate allele frequency from presence/absence

Calculate allele frequency from presence/absence

## Usage

``` r
calculate_af_presence_absence(allele_table)
```

## Arguments

- allele_table:

  Allele table with `target_name` and `variant`.

## Value

A tibble with `target_name`, `variant`, `allele_total`, `allele_count`,
and `freq`.
