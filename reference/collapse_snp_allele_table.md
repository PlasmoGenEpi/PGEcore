# Collapse overlapping-target SNP calls by summing or picking the best target

Collapse overlapping-target SNP calls by summing or picking the best
target

## Usage

``` r
collapse_snp_allele_table(
  allele_table_to_collapse,
  collapse_calls_by_summing = FALSE
)
```

## Arguments

- allele_table_to_collapse:

  SNP calls joined to allele counts.

- collapse_calls_by_summing:

  If `TRUE`, sum reads across targets.

## Value

Collapsed tibble.
