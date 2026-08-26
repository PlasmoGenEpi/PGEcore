# Collapse overlapping-target amino-acid calls

Collapse overlapping-target amino-acid calls

## Usage

``` r
collapse_aa_allele_table(
  allele_table_to_filter,
  collapse_calls_by_summing = FALSE
)
```

## Arguments

- allele_table_to_filter:

  Translated calls joined to allele counts.

- collapse_calls_by_summing:

  If `TRUE`, sum reads across targets.

## Value

Collapsed tibble.
