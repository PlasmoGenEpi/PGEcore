# Estimate allele frequency naively from AA or microhaplotype calls

Exactly one of `aa_calls` or `allele_table` must be provided. Frequency
is estimated either from within-sample read-count proportions
(`read_count_prop`) or from presence/absence (`presence_absence`).

## Usage

``` r
estimate_allele_frequency_naive(
  aa_calls = NULL,
  allele_table = NULL,
  method = "presence_absence",
  output = NULL
)
```

## Arguments

- aa_calls:

  Optional path to an AA calls TSV. See *Inputs*.

- allele_table:

  Optional path to an allele table TSV. See *Inputs*.

- method:

  Estimation method: `"presence_absence"` (default) or
  `"read_count_prop"`.

- output:

  Optional output TSV path. Default for the CLI is
  `allele_frequency.tsv`.

## Value

A tibble of estimated allele frequencies. For amino acid input the
`variant` column is a STAVE-style `gene_id:aa_position:aa` string. For
microhaplotype input columns include `target_name`, `seq`, and `freq`
(plus count columns for `presence_absence`).

## Details

### Inputs

- **`aa_calls`** (optional): AA calls (`specimen_name`, `gene_id`,
  `aa_position`, `aa`, `reads`). See
  [`vignette("input-formats", package = "PGEcore")`](https://plasmogenepi.github.io/PGEcore/articles/input-formats.md).

- **`allele_table`** (optional): Allele table (`specimen_name`,
  `target_name`, `seq`, `reads`). See the same vignette.

### Outputs

- **`output`** (optional): Allele-frequency TSV. For AA input, `variant`
  is a STAVE-style `gene_id:aa_position:aa` string plus `freq` (and
  count columns for `presence_absence`). For microhaplotype input:
  `target_name`, `seq`, `freq` (plus counts for `presence_absence`). If
  `NULL`, results are returned without writing a file.

### Running

    estimate_allele_frequency_naive(
      aa_calls = "aa_calls.tsv",
      output = "allele_frequency.tsv"
    )

    Rscript exec/estimate_allele_frequency_naive \
      --aa_calls aa_calls.tsv \
      --output allele_frequency.tsv

## See also

[`vignette("input-formats", package = "PGEcore")`](https://plasmogenepi.github.io/PGEcore/articles/input-formats.md)

## Examples

``` r
aa_path <- system.file(
  "extdata", "example_aa_calls.tsv",
  package = "PGEcore"
)
estimate_allele_frequency_naive(aa_calls = aa_path)
#> # A tibble: 10 × 4
#>    variant               allele_total allele_count  freq
#>    <chr>                        <int>        <int> <dbl>
#>  1 PF3D7_0417200.1:108:N            5            2 0.4  
#>  2 PF3D7_0417200.1:108:S            5            3 0.6  
#>  3 PF3D7_0417200.1:51:I             5            2 0.4  
#>  4 PF3D7_0417200.1:51:N             5            3 0.6  
#>  5 PF3D7_0417200.1:59:C             7            3 0.429
#>  6 PF3D7_0417200.1:59:R             7            4 0.571
#>  7 PF3D7_0810800.1:437:A            6            3 0.5  
#>  8 PF3D7_0810800.1:437:G            6            3 0.5  
#>  9 PF3D7_0810800.1:540:E            5            3 0.6  
#> 10 PF3D7_0810800.1:540:K            5            2 0.4  
```
