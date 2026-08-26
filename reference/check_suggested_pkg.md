# Check that a suggested package is installed

Used by wrappers that depend on optional R packages listed in Suggests.
Does not install the package.

## Usage

``` r
check_suggested_pkg(pkg, reason = NULL)
```

## Arguments

- pkg:

  Name of the package.

- reason:

  Optional short description of why it is needed.

## Value

Invisibly returns `TRUE` if the package is available.
