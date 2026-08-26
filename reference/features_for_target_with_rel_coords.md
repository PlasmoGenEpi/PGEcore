# Features intersecting a target, with coordinates relative to the target start

Features intersecting a target, with coordinates relative to the target
start

## Usage

``` r
features_for_target_with_rel_coords(lookup_row, feature_tab, intersect_col)
```

## Arguments

- lookup_row:

  One-row ref_bed tibble including the intersect column.

- feature_tab:

  Full feature table.

- intersect_col:

  Column of comma-separated feature row numbers.

## Value

Subset of `feature_tab` with `rel_start` and `rel_end`.
