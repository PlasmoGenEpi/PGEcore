# Mark which feature rows are fully contained in each panel interval

Mark which feature rows are fully contained in each panel interval

## Usage

``` r
add_intersected_features_to_ref_bed(ref_bed_tab, feature_tab, out_col)
```

## Arguments

- ref_bed_tab:

  Panel location table.

- feature_tab:

  Feature BED table (`#chrom`, `start`, `end`).

- out_col:

  Name of the column storing comma-separated feature row numbers.

## Value

`ref_bed_tab` with `out_col` added.
