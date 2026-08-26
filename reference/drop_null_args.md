# Drop NULL entries from a parsed optparse argument list

Used so
[`check_optparse_required_args()`](https://plasmogenepi.github.io/PGEcore/reference/check_optparse_required_args.md)
can treat unset options as missing.

## Usage

``` r
drop_null_args(arg)
```

## Arguments

- arg:

  Named list of parsed arguments.

## Value

`arg` with `NULL` values removed.
