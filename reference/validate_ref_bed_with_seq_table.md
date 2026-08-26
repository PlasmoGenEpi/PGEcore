# Validate a panel BED table that includes `ref_seq`

Validate a panel BED table that includes `ref_seq`

## Usage

``` r
validate_ref_bed_with_seq_table(ref_bed, data_name = "ref_bed")
```

## Arguments

- ref_bed:

  Data frame of panel locations with sequences.

- data_name:

  Label used in error messages.

## Value

Invisibly returns `TRUE` if valid.
