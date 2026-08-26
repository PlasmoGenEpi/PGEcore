# Summarize allele metrics by target name

Summarize allele metrics by target name

## Usage

``` r
summarize_allele_table(locus_data)
```

## Arguments

- locus_data:

  A data frame with columns `target_name` and `allele`.

## Value

A tibble with `target_name`, `total_allele_count`,
`unique_allele_count`, and `allele_singlets`.
