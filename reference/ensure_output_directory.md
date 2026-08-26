# Create an output directory, optionally replacing an existing one

Create an output directory, optionally replacing an existing one

## Usage

``` r
ensure_output_directory(output_directory, overwrite_dir = FALSE)
```

## Arguments

- output_directory:

  Directory path to create.

- overwrite_dir:

  If `TRUE`, delete `output_directory` first when it exists.

## Value

Invisibly returns `output_directory`.
