# Stop when validate::confront reports failing rules

Stop when validate::confront reports failing rules

## Usage

``` r
stop_on_validate_fails(df, rules, data_name)
```

## Arguments

- df:

  Data frame that was confronted.

- rules:

  A
  [`validate::validator()`](https://rdrr.io/pkg/validate/man/validator.html)
  object.

- data_name:

  Label used in the error message.

## Value

Invisibly returns `TRUE` if all rules pass.
