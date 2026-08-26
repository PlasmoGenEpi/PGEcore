# Check that an external executable is on PATH

External tools are not R package dependencies and are not installed by
PGEcore. Workflows must provide them in the process environment.

## Usage

``` r
check_external_tool(tool)
```

## Arguments

- tool:

  Executable name as found by
  [`Sys.which()`](https://rdrr.io/r/base/Sys.which.html).

## Value

Invisibly returns the absolute path to the executable.
