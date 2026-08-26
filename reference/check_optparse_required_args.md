# Stop if required optparse arguments are missing

Stop if required optparse arguments are missing

## Usage

``` r
check_optparse_required_args(arg, required_args)
```

## Arguments

- arg:

  Named list of parsed arguments (as from
  [`optparse::parse_args()`](https://trevorldavis.com/R/optparse/reference/parse_args.html)).

- required_args:

  Character vector of required argument names (without `--`).

## Value

Invisibly returns `TRUE` if all required arguments are present.
