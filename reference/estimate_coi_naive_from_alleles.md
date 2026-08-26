# Estimate COI from an allele-call data frame using naive methods

Estimate COI from an allele-call data frame using naive methods

## Usage

``` r
estimate_coi_naive_from_alleles(
  df_alleles,
  method = "integer_method",
  integer_threshold = 1,
  quantile_threshold = 0.05
)
```

## Arguments

- df_alleles:

  Data frame with columns `specimen_name`, `target_name`, `reads`, and
  `seq`.

- method:

  One of `integer_method` or `quantile_method`.

- integer_threshold:

  Index into the decreasing allele-count sequence (integer method only).

- quantile_threshold:

  Quantile in `[0, 1]` (quantile method only).

## Value

A tibble with columns `specimen_name` and `coi`.
