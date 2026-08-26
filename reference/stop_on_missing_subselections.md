# Stop when requested specimen/target names are missing from SNP data

Stop when requested specimen/target names are missing from SNP data

## Usage

``` r
stop_on_missing_subselections(
  snp_data,
  select_target_names,
  select_specimen_names,
  snp_table_fnp
)
```

## Arguments

- snp_data:

  SNP call tibble.

- select_target_names:

  Character vector of requested targets.

- select_specimen_names:

  Character vector of requested specimens.

- snp_table_fnp:

  Path label used in messages.

## Value

Invisibly `TRUE` if all requested names are present.
