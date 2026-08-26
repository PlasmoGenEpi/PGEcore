# Estimate allele prevalence naively from AA or microhaplotype calls

Exactly one of `aa_calls` or `allele_table` must be provided. Prevalence
is the fraction of specimens carrying each allele at a locus.

## Usage

``` r
estimate_allele_prevalence_naive(
  aa_calls = NULL,
  allele_table = NULL,
  output = NULL
)
```

## Arguments

- aa_calls:

  Optional path to an AA calls TSV. See *Inputs*.

- allele_table:

  Optional path to an allele table TSV. See *Inputs*.

- output:

  Optional output TSV path. Default for the CLI is `prevalence.tsv`.

## Value

A tibble of allele prevalences. For amino acid input: `variant`, `prev`,
`sample_count`, `sample_total`. For microhaplotype input: `target_name`,
`seq`, `prev`, `sample_count`, `sample_total`.

## Details

### Inputs

- **`aa_calls`** (optional): AA calls (`specimen_name`, `gene_id`,
  `aa_position`, `aa`). See
  [`vignette("input-formats", package = "PGEcore")`](https://plasmogenepi.github.io/PGEcore/articles/input-formats.md).

- **`allele_table`** (optional): Allele table (`specimen_name`,
  `target_name`, `seq`). See the same vignette.

### Outputs

- **`output`** (optional): Prevalence TSV. For AA input: `variant`,
  `prev`, `sample_count`, `sample_total`. For microhaplotype input:
  `target_name`, `seq`, `prev`, `sample_count`, `sample_total`. If
  `NULL`, results are returned without writing a file.

### Running

    estimate_allele_prevalence_naive(
      aa_calls = "aa_calls.tsv",
      output = "prevalence.tsv"
    )

    Rscript exec/estimate_allele_prevalence_naive \
      --aa_calls aa_calls.tsv \
      --output prevalence.tsv

## See also

[`vignette("input-formats", package = "PGEcore")`](https://plasmogenepi.github.io/PGEcore/articles/input-formats.md)

## Examples

``` r
aa_path <- system.file(
  "extdata", "example_aa_calls.tsv",
  package = "PGEcore"
)
estimate_allele_prevalence_naive(aa_calls = aa_path)
#> # A tibble: 10 × 4
#>    variant                prev sample_count sample_total
#>    <chr>                 <dbl>        <int>        <int>
#>  1 PF3D7_0417200.1:108:N  0.5             2            4
#>  2 PF3D7_0417200.1:108:S  0.75            3            4
#>  3 PF3D7_0417200.1:51:I   0.5             2            4
#>  4 PF3D7_0417200.1:51:N   0.75            3            4
#>  5 PF3D7_0417200.1:59:C   0.75            3            4
#>  6 PF3D7_0417200.1:59:R   1               4            4
#>  7 PF3D7_0810800.1:437:A  0.75            3            4
#>  8 PF3D7_0810800.1:437:G  0.75            3            4
#>  9 PF3D7_0810800.1:540:E  0.75            3            4
#> 10 PF3D7_0810800.1:540:K  0.5             2            4
```
