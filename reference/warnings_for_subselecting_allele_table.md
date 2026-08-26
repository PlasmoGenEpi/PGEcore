# Warnings when requested specimen/target names are missing from an allele table

Warnings when requested specimen/target names are missing from an allele
table

## Usage

``` r
warnings_for_subselecting_allele_table(
  allele_data,
  select_target_names,
  select_specimen_names,
  allele_table_fnp
)
```

## Arguments

- allele_data:

  Allele table.

- select_target_names:

  Requested targets (empty = none).

- select_specimen_names:

  Requested specimens (empty = none).

- allele_table_fnp:

  Path label used in messages.

## Value

Character vector of warning strings (possibly empty).
