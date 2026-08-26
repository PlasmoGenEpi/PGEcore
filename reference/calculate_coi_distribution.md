# Calculate the distribution of COI values across specimens

Calculate the distribution of COI values across specimens

## Usage

``` r
calculate_coi_distribution(coi_table)
```

## Arguments

- coi_table:

  A data frame with a numeric `coi` column.

## Value

A tibble with columns `coi`, `n`, and `proportion` for each integer COI
from 1 to `max(coi)`.
