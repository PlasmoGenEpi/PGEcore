# Collect validate::confront failures as a warning string (or NULL)

Matches legacy scripts that recorded validation problems in `warns` and
continued rather than stopping.

## Usage

``` r
warn_on_validate_fails(df, rules, data_name)
```

## Arguments

- df:

  Data frame that was confronted.

- rules:

  A
  [`validate::validator()`](https://rdrr.io/pkg/validate/man/validator.html)
  object.

- data_name:

  Label used in the message.

## Value

Character warning message, or `NULL` if all rules pass.
