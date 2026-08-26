# Estimate multilocus prevalence and frequency with naive phasing

For each loci group, specimens that have calls at every locus are
retained. Haplotypes are inferred when exactly one locus is
heterozygous, or when every locus has a single allele above
`wsaf_cut_off` (including monoclonal samples). Prevalence and frequency
are then estimated with `wsaf_prop` or `presence_absence`.

## Usage

``` r
multilocus_prevfreq_naive(
  aa_calls,
  loci_groups,
  output = NULL,
  single_locus_output = NULL,
  method = "wsaf_prop",
  wsaf_cut_off = 0.7
)
```

## Arguments

- aa_calls:

  Path to an AA calls TSV. See *Inputs*.

- loci_groups:

  Path to a loci groups TSV. See *Inputs*.

- output:

  Optional path for the multilocus prev/freq TSV.

- single_locus_output:

  Optional path for single-locus prev/freq recalculated from the
  inferred multilocus calls.

- method:

  `"wsaf_prop"` (default) or `"presence_absence"`.

- wsaf_cut_off:

  WSAF threshold used to infer a dominant haplotype when more than one
  locus is heterozygous. Default: `0.70`.

## Value

A tibble of multilocus prevalence and frequency estimates.

## Details

### Inputs

- **`aa_calls`**: AA calls (`specimen_name`, `gene`, `gene_id`,
  `aa_position`, `reads`, `aa`). See
  [`vignette("input-formats", package = "PGEcore")`](https://plasmogenepi.github.io/PGEcore/articles/input-formats.md).

- **`loci_groups`**: Loci groups (`group_id`, `gene_id`, `aa_position`).
  See the same vignette.

### Outputs

- **`output`** (optional): Multilocus prev/freq TSV (includes
  `group_id`, `variant`, `prev`, `freq`). If `NULL`, results are
  returned without writing.

- **`single_locus_output`** (optional): Single-locus prev/freq
  recalculated from the inferred multilocus calls.

### Running

    multilocus_prevfreq_naive(
      aa_calls = "aa_calls.tsv",
      loci_groups = "loci_groups.tsv",
      output = "multilocus_prevfreq.tsv"
    )

    Rscript exec/multilocus_prevfreq_naive \
      --aa_calls aa_calls.tsv \
      --loci_groups loci_groups.tsv \
      --output multilocus_prevfreq.tsv

## See also

[`vignette("input-formats", package = "PGEcore")`](https://plasmogenepi.github.io/PGEcore/articles/input-formats.md)

## Examples

``` r
aa_path <- system.file(
  "extdata", "example2_aa_calls.tsv",
  package = "PGEcore"
)
groups_path <- system.file(
  "extdata", "example_loci_groups.tsv",
  package = "PGEcore"
)
multilocus_prevfreq_naive(aa_path, groups_path)
#> # A tibble: 12 × 6
#>    variant                      sample_total sample_count   prev   freq group_id
#>    <chr>                               <int>        <int>  <dbl>  <dbl> <chr>   
#>  1 PF3D7_0417200.1:51_59_108:I…           25           18 0.72   0.590  pfdhfr  
#>  2 PF3D7_0417200.1:51_59_108:N…           25           12 0.48   0.410  pfdhfr  
#>  3 PF3D7_0417200.1:51_59_108:I…           24            4 0.167  0.170  pfdhfr_…
#>  4 PF3D7_0417200.1:51_59_108:I…           24            3 0.125  0.0798 pfdhfr_…
#>  5 PF3D7_0417200.1:51_59_108:I…           24            4 0.167  0.163  pfdhfr_…
#>  6 PF3D7_0417200.1:51_59_108:I…           24            4 0.167  0.170  pfdhfr_…
#>  7 PF3D7_0417200.1:51_59_108:N…           24            8 0.333  0.332  pfdhfr_…
#>  8 PF3D7_0417200.1:51_59_108:N…           24            2 0.0833 0.0851 pfdhfr_…
#>  9 PF3D7_0810800.1:437_540:A_K            24           12 0.5    0.502  pfdhps  
#> 10 PF3D7_0810800.1:437_540:G_E            24            3 0.125  0.0797 pfdhps  
#> 11 PF3D7_0810800.1:437_540:G_K            24            7 0.292  0.249  pfdhps  
#> 12 PF3D7_0810800.1:437_540:G_N            24            4 0.167  0.170  pfdhps  
```
