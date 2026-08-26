# Read and validate a locus table for per-locus popgen summaries

Read and validate a locus table for per-locus popgen summaries

## Usage

``` r
read_popgen_locus_data(
  input_path,
  specimen_name_col = "specimen_name",
  target_name_col = "target_name",
  target_value_col = "seq"
)
```

## Arguments

- input_path:

  Path to a TSV allele table.

- specimen_name_col:

  Specimen ID column name.

- target_name_col:

  Target/locus column name.

- target_value_col:

  Allele/sequence column name.

## Value

A tibble with `specimen_name`, `target_name`, and `target_value`.
