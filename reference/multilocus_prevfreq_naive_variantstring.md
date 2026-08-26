# Estimate multilocus prevalence and frequency with variantstring

Converts amino acid calls to [variant
strings](https://github.com/mrc-ide/variantstring), extracts unambiguous
component genotypes, and estimates prevalence and frequency of each
component against the full dataset.

## Usage

``` r
multilocus_prevfreq_naive_variantstring(aa_calls, loci_groups, output = NULL)
```

## Arguments

- aa_calls:

  Path to an AA calls TSV. See *Inputs*.

- loci_groups:

  Path to a loci groups TSV. See *Inputs*.

- output:

  Optional path for the prev/freq TSV.

## Value

A tibble with columns `group_id`, `variant`, `prev`, `freq`, and
`sample_total`.

## Details

### Inputs

- **`aa_calls`**: AA calls (`specimen_name`, `gene_id`, `aa_position`,
  `reads`, `aa`). See
  [`vignette("input-formats", package = "PGEcore")`](https://plasmogenepi.github.io/PGEcore/articles/input-formats.md).

- **`loci_groups`**: Loci groups (`group_id`, `gene_id`, `aa_position`).
  See the same vignette.

### Outputs

- **`output`** (optional): Prev/freq TSV with columns `group_id`,
  `variant`, `prev`, `freq`, and `sample_total`. If `NULL`, results are
  returned without writing.

### Running

    multilocus_prevfreq_naive_variantstring(
      aa_calls = "aa_calls.tsv",
      loci_groups = "loci_groups.tsv",
      output = "multilocus_prevfreq.tsv"
    )

    Rscript exec/multilocus_prevfreq_naive_variantstring \
      --aa_calls aa_calls.tsv \
      --loci_groups loci_groups.tsv \
      --output multilocus_prevfreq.tsv

Requires the optional **variantstring** package (Suggests). It is not
installed automatically with PGEcore.

## See also

[`vignette("input-formats", package = "PGEcore")`](https://plasmogenepi.github.io/PGEcore/articles/input-formats.md)

## Examples

``` r
aa_path <- system.file(
  "extdata", "example_aa_calls.tsv",
  package = "PGEcore"
)
groups_path <- system.file(
  "extdata", "example_loci_groups.tsv",
  package = "PGEcore"
)
multilocus_prevfreq_naive_variantstring(aa_path, groups_path)
#> # A tibble: 12 × 5
#>    group_id      variant                                prev   freq sample_total
#>    <chr>         <chr>                                 <dbl>  <dbl>        <int>
#>  1 pfdhfr_pfdhps PF3D7_0417200.1:51_59_108:N_R_N;PF3D… 0.333 0.333             3
#>  2 pfdhfr_pfdhps PF3D7_0417200.1:51_59_108:N_C_S;PF3D… 0.333 0.267             3
#>  3 pfdhfr_pfdhps PF3D7_0417200.1:51_59_108:N_R_S;PF3D… 0.333 0.0667            3
#>  4 pfdhfr        PF3D7_0417200.1:51_59_108:N_R_N       0.333 0.333             3
#>  5 pfdhfr        PF3D7_0417200.1:51_59_108:N_C_S       0.333 0.267             3
#>  6 pfdhfr        PF3D7_0417200.1:51_59_108:N_R_S       0.333 0.0667            3
#>  7 pfdhfr        PF3D7_0417200.1:51_59_108:I_C_S       0.333 0.133             3
#>  8 pfdhfr        PF3D7_0417200.1:51_59_108:I_R_S       0.333 0.2               3
#>  9 pfdhps        PF3D7_0810800.1:437_540:A_E           0.333 0.333             3
#> 10 pfdhps        PF3D7_0810800.1:437_540:G_E           0.333 0.333             3
#> 11 pfdhps        PF3D7_0810800.1:437_540:A_K           0.333 0.0667            3
#> 12 pfdhps        PF3D7_0810800.1:437_540:G_K           0.333 0.267             3
```
