# Count specimens by complexity of infection (COI)

Reads per-specimen COI values, rounds them to integers, and returns the
count and proportion of specimens at each COI level.

## Usage

``` r
count_samples_by_coi(coi_table, output = NULL)
```

## Arguments

- coi_table:

  Path to a COI table TSV, or a data frame with the same columns. See
  *Inputs*.

- output:

  Optional output TSV path. Default for the CLI is
  `coi_distribution.tsv`.

## Value

A tibble with columns `coi`, `n`, and `proportion`.

## Details

### Inputs

- **`coi_table`**: COI table (`specimen_name`, `coi`), as a file path or
  data frame. See
  [`vignette("input-formats", package = "PGEcore")`](https://plasmogenepi.github.io/PGEcore/articles/input-formats.md).

### Outputs

- **`output`** (optional): TSV with columns `coi`, `n`, and
  `proportion`. If `NULL`, results are returned without writing a file.

### Running

    count_samples_by_coi(
      coi_table = "coi_table.tsv",
      output = "coi_distribution.tsv"
    )

    Rscript exec/count_samples_by_coi \
      --coi_table coi_table.tsv \
      --output coi_distribution.tsv

## See also

[`vignette("input-formats", package = "PGEcore")`](https://plasmogenepi.github.io/PGEcore/articles/input-formats.md)

## Examples

``` r
coi_path <- system.file("extdata", "example_coi_table.tsv", package = "PGEcore")
count_samples_by_coi(coi_path)
#> # A tibble: 3 × 3
#>     coi     n proportion
#>   <dbl> <int>      <dbl>
#> 1     1     0        0  
#> 2     2     1        0.2
#> 3     3     4        0.8
```
