# Filter an SNP table by optional target and specimen selections

Filter an SNP table by optional target and specimen selections

## Usage

``` r
filter_snp_table_for_optional_subselecting(
  snp_data,
  select_target_names = character(0),
  select_specimen_names = character(0)
)
```

## Arguments

- snp_data:

  SNP call tibble.

- select_target_names:

  Character vector of target names (empty = no filter).

- select_specimen_names:

  Character vector of specimen names (empty = no filter).

## Value

Filtered `snp_data`.
