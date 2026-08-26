# Read amino acid calls for variantstring multilocus prev/freq

Read amino acid calls for variantstring multilocus prev/freq

## Usage

``` r
read_mlp_vs_aa_table(aa_calls)
```

## Arguments

- aa_calls:

  Path to a TSV of amino acid calls.

## Value

A tibble with `specimen_name`, `gene`, `pos`, `reads`, `aa`, `n_aa`.
