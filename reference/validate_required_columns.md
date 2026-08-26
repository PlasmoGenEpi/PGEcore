# Validate that a data frame has required columns and is non-empty

Validate that a data frame has required columns and is non-empty

## Usage

``` r
validate_required_columns(data, required_cols, data_name)
```

## Arguments

- data:

  Data frame to validate.

- required_cols:

  Character vector of required column names.

- data_name:

  Label used in error messages.

## Value

Invisibly returns `TRUE` if valid.
