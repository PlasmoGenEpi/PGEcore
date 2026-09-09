# Estimate single-locus allele frequencies with the Incomplete Data Model

Estimates single-locus allele frequencies with the Incomplete Data Model
(or original model). Provide exactly one of `allele_table` or
`aa_calls`. Vendored MLE code is from Hashemi & Schneider (2024).
Requires **Rmpfr** and **openxlsx** (Suggests).

## Usage

``` r
IDM_wrapper(
  allele_table = "",
  aa_calls = "",
  slaf_output,
  model = "IDM",
  lambda_initial = 1,
  eps_initial = 0.1
)
```

## Arguments

- allele_table:

  Path to allele table TSV, or `""` / `NULL` if using amino-acid calls.
  See *Inputs*.

- aa_calls:

  Path to amino-acid calls TSV, or `""` / `NULL` if using an allele
  table. See *Inputs*.

- slaf_output:

  Output TSV path. See *Outputs*.

- model:

  `"IDM"` (incomplete-data model), `"OM"` (original model), or
  `"IDM_OM"` (IDM, falling back to the OM at any locus the IDM leaves
  unsolved), `"OM_CC"` (OM with a continuity correction on lineage
  prevalence), or `"IDM_OM_CC"` (IDM, falling back to `"OM_CC"`).

- lambda_initial:

  Initial lambda for the numerical iteration.

- eps_initial:

  Initial epsilon for the numerical iteration.

## Value

The SLAF tibble (invisibly after writing `slaf_output`).

## Details

### Inputs

- **`allele_table`**: Allele table TSV, or `""` / `NULL` if using
  `aa_calls`. See
  [`vignette("input-formats", package = "PGEcore")`](https://plasmogenepi.github.io/PGEcore/articles/input-formats.md).

- **`aa_calls`**: Amino-acid calls TSV, or `""` / `NULL` if using
  `allele_table`.

### Outputs

- **`slaf_output`**: Single-locus allele frequencies. Allele-table input
  is written as `target_name`, `seq`, `freq`; AA-call input as
  `variant`, `freq`.

### Running

    IDM_wrapper(
      allele_table = "allele_table.tsv",
      slaf_output = "slaf.tsv"
    )

    Rscript exec/IDM_wrapper \
      --allele_table allele_table.tsv \
      --slaf_output slaf.tsv

Requires **Rmpfr** and **openxlsx** (Suggests).

## See also

[`vignette("input-formats", package = "PGEcore")`](https://plasmogenepi.github.io/PGEcore/articles/input-formats.md)
