# Per-locus nucleotide diversity, segregating sites, and Tajima's D

Groups allele sequences by locus and computes population-genetic
summaries.

## Usage

``` r
per_locus_popgen_summary(
  allele_table,
  specimen_name_col = "specimen_name",
  target_name_col = "target_name",
  target_value_col = "seq",
  output = "per_locus_popgen_summary.tsv",
  msa_method = "Muscle"
)
```

## Arguments

- allele_table:

  Path to an allele table TSV. See *Inputs*.

- specimen_name_col:

  Specimen ID column name.

- target_name_col:

  Locus column name.

- target_value_col:

  Allele/sequence column name.

- output:

  Optional output TSV path. Defaults to
  `"per_locus_popgen_summary.tsv"`.

- msa_method:

  Alignment method: `"Muscle"` (default), `"ClustalW"`, or
  `"ClustalOmega"`.

## Value

A tibble of per-locus statistics with lower-case column names.

## Details

### Inputs

- **`allele_table`**: Allele table (default columns `specimen_name`,
  `target_name`, `seq`). See
  [`vignette("input-formats", package = "PGEcore")`](https://plasmogenepi.github.io/PGEcore/articles/input-formats.md).

### Outputs

- **`output`** (optional): Per-locus stats TSV with lower-case column
  names (e.g. `target_name`, `nucleotide_diversity`,
  `segregating_sites`, `tajima_d`, …). Defaults to
  `"per_locus_popgen_summary.tsv"`. If `NULL`, results are returned
  without writing.

### Running

    per_locus_popgen_summary(
      allele_table = "allele_table.tsv",
      output = "per_locus_popgen_summary.tsv"
    )

    Rscript exec/per_locus_popgen_summary \
      --allele_table allele_table.tsv \
      --output per_locus_popgen_summary.tsv

Requires **ape**, **msa**, and **pegas** (Suggests). **msa** calls an
external aligner; install the binary for `msa_method` on `PATH`:

- `"Muscle"` — `muscle`

- `"ClustalW"` — `clustalw`

- `"ClustalOmega"` — `clustalo`

## See also

[`vignette("input-formats", package = "PGEcore")`](https://plasmogenepi.github.io/PGEcore/articles/input-formats.md)
