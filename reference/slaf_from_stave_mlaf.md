# Convert STAVE multi-locus allele frequencies to single-locus frequencies

Expands STAVE `variant` strings with **variantstring**, aggregates
frequencies per amino acid allele, and emits STAVE-style single-locus
`variant` identifiers via
[`convert_single_locus_table_to_stave()`](https://plasmogenepi.github.io/PGEcore/reference/convert_single_locus_table_to_stave.md).

## Usage

``` r
slaf_from_stave_mlaf(mlaf, output = NULL)
```

## Arguments

- mlaf:

  Path to an MLAF TSV, or a data frame with the same columns. See
  *Inputs*.

- output:

  Optional output TSV path. Default for the CLI is
  `single_locus_allele_frequencies.tsv`.

## Value

A tibble with columns `variant` and `freq` (legacy STAVE conversion
drops `group_id`, matching the original script).

## Details

### Inputs

- **`mlaf`**: Multilocus allele frequency table (`group_id`, `variant`,
  `freq`), as a file path or data frame. See
  [`vignette("input-formats", package = "PGEcore")`](https://plasmogenepi.github.io/PGEcore/articles/input-formats.md).

### Outputs

- **`output`** (optional): Single-locus allele frequency TSV with
  columns `variant` and `freq`. If `NULL`, results are returned without
  writing. CLI default is `single_locus_allele_frequencies.tsv`.

### Running

    slaf_from_stave_mlaf(
      mlaf = "mlaf.tsv",
      output = "single_locus_allele_frequencies.tsv"
    )

    Rscript exec/slaf_from_stave_mlaf \
      --mlaf mlaf.tsv \
      --output single_locus_allele_frequencies.tsv

Requires the optional **variantstring** package (Suggests). It is not
installed automatically with PGEcore.

## See also

[`vignette("input-formats", package = "PGEcore")`](https://plasmogenepi.github.io/PGEcore/articles/input-formats.md)
