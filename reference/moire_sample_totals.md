# Number of specimens with a call at each locus

Number of specimens with a call at each locus

## Usage

``` r
moire_sample_totals(moire_data)
```

## Arguments

- moire_data:

  Long-form table with `sample_id`, `locus` and `allele`.

## Value

A tibble with `target_name` and `sample_total`, the number of distinct
specimens observed at that locus.
