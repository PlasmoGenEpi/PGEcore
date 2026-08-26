# Summarize alleles per locus from an allele table

For each `target_name`, computes total allele count, unique allele
count, and the number of alleles that appear only once (singlets).

## Usage

``` r
allele_per_locus_summary(allele_table, output = "allele_summary_by_target.tsv")
```

## Arguments

- allele_table:

  Path to an allele table TSV. See *Inputs*.

- output:

  Optional output TSV path. Defaults to
  `"allele_summary_by_target.tsv"`.

## Value

A tibble with columns `target_name`, `total_allele_count`,
`unique_allele_count`, and `allele_singlets`.

## Details

### Inputs

- **`allele_table`**: Allele table (`specimen_name`, `target_name`,
  `seq`). See
  [`vignette("input-formats", package = "PGEcore")`](https://plasmogenepi.github.io/PGEcore/articles/input-formats.md).

### Outputs

- **`output`** (optional): TSV with columns `target_name`,
  `total_allele_count`, `unique_allele_count`, and `allele_singlets`.
  Defaults to `"allele_summary_by_target.tsv"`. If `NULL`, results are
  returned without writing a file.

### Running

    allele_per_locus_summary(
      allele_table = "allele_table.tsv",
      output = "allele_summary_by_target.tsv"
    )

    Rscript exec/allele_per_locus_summary \
      --allele_table allele_table.tsv \
      --output allele_summary_by_target.tsv

## See also

[`vignette("input-formats", package = "PGEcore")`](https://plasmogenepi.github.io/PGEcore/articles/input-formats.md)

## Examples

``` r
allele_path <- system.file(
  "extdata", "example_allele_table.tsv",
  package = "PGEcore"
)
allele_per_locus_summary(allele_path, output = NULL)
#> [1] "Reading input data"
#> [1] "Validating input format"
#> [1] "Confronting input data with validation rules"
#> [1] "Returning Locus data"
#> # A tibble: 5 × 4
#>   target_name             total_allele_count unique_allele_count allele_singlets
#>   <chr>                                <int>               <int>           <int>
#> 1 Pf3D7_01_v3-145388-145…                 11                   3               0
#> 2 Pf3D7_01_v3-162867-163…                 12                   3               0
#> 3 Pf3D7_01_v3-181512-181…                 11                   3               0
#> 4 Pf3D7_01_v3-194742-194…                 10                   3               0
#> 5 Pf3D7_01_v3-455794-456…                  6                   2               1
```
