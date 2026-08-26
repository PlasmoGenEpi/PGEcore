# Stop if an output path exists and overwrite is FALSE

Stop if an output path exists and overwrite is FALSE

## Usage

``` r
stop_if_output_exists(path, overwrite = FALSE)
```

## Arguments

- path:

  Output file path (ignored when `NULL`).

- overwrite:

  Whether overwriting is allowed.

## Value

Invisibly returns `TRUE` if writing may proceed.
