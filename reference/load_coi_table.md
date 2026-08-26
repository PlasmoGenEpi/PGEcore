# Load COI calls from a TSV file

Reads `specimen_name` and `coi`, and rounds `coi` to the nearest
integer.

## Usage

``` r
load_coi_table(path)
```

## Arguments

- path:

  Path to a TSV with columns `specimen_name` and `coi`.

## Value

A tibble with columns `specimen_name` and integer-rounded `coi`.
