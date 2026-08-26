# Population-genetic statistics grouped by `target_name`

Population-genetic statistics grouped by `target_name`

## Usage

``` r
calculate_stats_by_target_name(locus_data, msa_method = "Muscle")
```

## Arguments

- locus_data:

  Tibble with `target_name` and `target_value`.

- msa_method:

  Alignment method for
  [`calculate_popgen_stats()`](https://plasmogenepi.github.io/PGEcore/reference/calculate_popgen_stats.md).

## Value

A tibble of per-locus statistics.
