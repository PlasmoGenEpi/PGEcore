# Add `covered_by_target` listing panel targets that fully contain each feature

Add `covered_by_target` listing panel targets that fully contain each
feature

## Usage

``` r
add_covered_by_target_to_features(feature_tab, ref_bed_tab)
```

## Arguments

- feature_tab:

  Feature BED table.

- ref_bed_tab:

  Panel location table.

## Value

`feature_tab` with `covered_by_target` (`"uncovered"` if none).
