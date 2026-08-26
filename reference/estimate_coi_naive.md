# Estimate COI using naive allele-count methods

For every specimen, counts distinct alleles at each locus and sorts
those counts in decreasing order. With `method = "integer_method"`, the
`integer_threshold`-th value is the COI estimate. With
`method = "quantile_method"`, the value at `quantile_threshold` is used
instead (scaling naturally with the number of observed allele rows).

## Usage

``` r
estimate_coi_naive(
  allele_table,
  output = NULL,
  method = "integer_method",
  integer_threshold = 1,
  quantile_threshold = 0.05
)
```

## Arguments

- allele_table:

  Path to an allele table TSV. See *Inputs*.

- output:

  Optional output TSV path.

- method:

  One of `"integer_method"` or `"quantile_method"`. Default:
  `"integer_method"`.

- integer_threshold:

  Positive integer index into the ordered allele counts (integer method
  only). Default: `1`.

- quantile_threshold:

  Quantile in `[0, 1]` (quantile method only). Values near zero yield
  higher COI estimates. Default: `0.05`.

## Value

A tibble with columns `specimen_name` and `coi`.

## Details

### Inputs

- **`allele_table`**: Allele table (`specimen_name`, `target_name`,
  `seq`, `reads`). See
  [`vignette("input-formats", package = "PGEcore")`](https://plasmogenepi.github.io/PGEcore/articles/input-formats.md).

### Outputs

- **`output`** (optional): COI table TSV with columns `specimen_name`
  and `coi`. If `NULL`, results are returned without writing a file.

### Running

    estimate_coi_naive(
      allele_table = "allele_table.tsv",
      output = "coi_table.tsv"
    )

    Rscript exec/estimate_coi_naive \
      --allele_table allele_table.tsv \
      --output coi_table.tsv

## See also

[`vignette("input-formats", package = "PGEcore")`](https://plasmogenepi.github.io/PGEcore/articles/input-formats.md)

## Examples

``` r
allele_path <- system.file(
  "extdata", "example_allele_table.tsv",
  package = "PGEcore"
)
estimate_coi_naive(allele_path, method = "integer_method")
#> # A tibble: 5 × 2
#>   specimen_name                                                              coi
#>   <chr>                                                                    <int>
#> 1 PARAV3-ENV-MH04-7S1-7C1-1000-parasitedensity-sampleDB-4064911117_S43_L0…     3
#> 2 PARAV3-ENV-MH04-DS2-DC11-1000-parasitedensity-sampleDB-4064911661_S35_L…     2
#> 3 PARAV3-ENV-MH04-DS2-DC11-10000-parasitedensity-sampleDB-4064911565_S30_…     3
#> 4 PARAV3-ENV-MH04-DS2-DC3-1000-parasitedensity-sampleDB-4064911200_S31_L0…     3
#> 5 PARAV3-ENV-MH04-DS4-DC2-1000-parasitedensity-sampleDB-4064911921_S14_L0…     3
```
