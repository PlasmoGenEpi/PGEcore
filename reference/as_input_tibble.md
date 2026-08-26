# Coerce ref_bed / allele_table / feature table inputs to tibbles

Coerce ref_bed / allele_table / feature table inputs to tibbles

## Usage

``` r
as_input_tibble(x, reader, what)
```

## Arguments

- x:

  A path or data frame.

- reader:

  Function used when `x` is a path.

- what:

  Label for error messages.

## Value

A tibble.
