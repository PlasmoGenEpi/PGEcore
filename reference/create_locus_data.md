# Create locus data from an allele table file

Reads a TSV with columns `specimen_name`, `target_name`, and `seq`,
renames them to `sample_id`, `target_name`, and `allele`, and validates
character non-missing values with the **validate** package.

## Usage

``` r
create_locus_data(input_path)
```

## Arguments

- input_path:

  Path to the allele table TSV.

## Value

A data frame with columns `sample_id`, `target_name`, and `allele`.
