# Stop if a character vector contains duplicated values

Stop if a character vector contains duplicated values

## Usage

``` r
stop_on_duplicate_names(values, source, field = "target_name")
```

## Arguments

- values:

  Values to check.

- source:

  Label for the input file or table.

- field:

  Field name used in the error message.

## Value

Invisibly returns `TRUE` if all values are unique.
