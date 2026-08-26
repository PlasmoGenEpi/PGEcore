# Convert a single-locus table to STAVE-style variant identifiers

Builds a `variant` column by concatenating `gene_id`, `aa_position`, and
`aa` as `gene_id:aa_position:aa`, then returns `variant` plus any
requested additional columns.

## Usage

``` r
convert_single_locus_table_to_stave(df, additional_columns = NULL)
```

## Arguments

- df:

  A data frame with columns `gene_id`, `aa_position`, and `aa`, plus any
  columns named in `additional_columns`.

- additional_columns:

  Optional character vector of extra columns to keep.

## Value

A data frame with a `variant` column and any `additional_columns`.

## Examples

``` r
df <- data.frame(
  gene_id = c("PF3D7_0417200.1", "PF3D7_0417200.1"),
  aa_position = c(51, 59),
  aa = c("I", "R"),
  prev = c(0.5, 1)
)
convert_single_locus_table_to_stave(df, additional_columns = "prev")
#>                variant prev
#> 1 PF3D7_0417200.1:51:I  0.5
#> 2 PF3D7_0417200.1:59:R  1.0
```
