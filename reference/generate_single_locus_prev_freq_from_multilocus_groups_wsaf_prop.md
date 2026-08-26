# Recalculate single-locus prev/freq from multilocus calls (WSAF-weighted)

Recalculate single-locus prev/freq from multilocus calls (WSAF-weighted)

## Usage

``` r
generate_single_locus_prev_freq_from_multilocus_groups_wsaf_prop(
  multilocus_calls
)
```

## Arguments

- multilocus_calls:

  Table with `group_id`, `specimen_name`, `variant`, `wsaf`.

## Value

Summarized single-locus prev/freq per group.
