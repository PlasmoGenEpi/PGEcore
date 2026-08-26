# Filter amino acid calls to biallelic loci

Keeps loci (`gene_id`, `aa_position`, `ref_aa`) with at most two
distinct `aa` alleles. Optionally writes non-biallelic loci to a second
file.

## Usage

``` r
filter_biallelic_calls(
  aa_calls,
  output = NULL,
  nonbiallelic_output = NULL,
  overwrite = FALSE
)
```

## Arguments

- aa_calls:

  Path to an AA calls TSV, or a data frame with the same columns. See
  *Inputs*.

- output:

  Optional path for the biallelic output TSV.

- nonbiallelic_output:

  Optional path for loci with more than two alleles.

- overwrite:

  If `FALSE` (default), refuse to overwrite existing outputs.

## Value

A list with tibbles `biallelic` and `nonbiallelic`. Each includes an
`allele_calls` column with the distinct allele count per locus.

## Details

### Inputs

- **`aa_calls`**: AA calls (`gene_id`, `aa_position`, `ref_aa`, `aa`),
  as a file path or data frame. See
  [`vignette("input-formats", package = "PGEcore")`](https://plasmogenepi.github.io/PGEcore/articles/input-formats.md).

### Outputs

- **`output`** (optional): Biallelic AA calls TSV (includes
  `allele_calls`).

- **`nonbiallelic_output`** (optional): Non-biallelic loci TSV (includes
  `allele_calls`).

### Running

    filter_biallelic_calls(
      aa_calls = "aa_calls.tsv",
      output = "biallelic_aa_calls.tsv"
    )

    Rscript exec/filter_biallelic_calls \
      --aa_calls aa_calls.tsv \
      --output biallelic_aa_calls.tsv

## See also

[`vignette("input-formats", package = "PGEcore")`](https://plasmogenepi.github.io/PGEcore/articles/input-formats.md)

## Examples

``` r
path <- system.file("extdata", "example_aa_calls.tsv", package = "PGEcore")
filter_biallelic_calls(path)
#> $biallelic
#> # A tibble: 28 × 10
#>    specimen_name target_name    reads gene   aa_locus gene_id aa_position ref_aa
#>    <chr>         <chr>          <dbl> <chr>  <chr>    <chr>         <dbl> <chr> 
#>  1 specimen1     pfdhfr_1_150       5 dhfr-… PF3D7_0… PF3D7_…          51 N     
#>  2 specimen1     pfdhfr_1_150       5 dhfr-… PF3D7_0… PF3D7_…          59 C     
#>  3 specimen1     pfdhfr_1_150       5 dhfr-… PF3D7_0… PF3D7_…         108 S     
#>  4 specimen1     pfdhps_400_550     5 dhps   PF3D7_0… PF3D7_…         437 A     
#>  5 specimen1     pfdhps_400_550     5 dhps   PF3D7_0… PF3D7_…         540 K     
#>  6 specimen2     pfdhfr_1_150       5 dhfr-… PF3D7_0… PF3D7_…          51 N     
#>  7 specimen2     pfdhfr_1_150       4 dhfr-… PF3D7_0… PF3D7_…          59 C     
#>  8 specimen2     pfdhfr_1_150       1 dhfr-… PF3D7_0… PF3D7_…          59 C     
#>  9 specimen2     pfdhfr_1_150       5 dhfr-… PF3D7_0… PF3D7_…         108 S     
#> 10 specimen2     pfdhps_400_550     5 dhps   PF3D7_0… PF3D7_…         437 A     
#> # ℹ 18 more rows
#> # ℹ 2 more variables: aa <chr>, allele_calls <int>
#> 
#> $nonbiallelic
#> # A tibble: 0 × 10
#> # ℹ 10 variables: specimen_name <chr>, target_name <chr>, reads <dbl>,
#> #   gene <chr>, aa_locus <chr>, gene_id <chr>, aa_position <dbl>, ref_aa <chr>,
#> #   aa <chr>, allele_calls <int>
#> 
```
