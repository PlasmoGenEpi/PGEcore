# Calculate single-locus allele frequencies from microhaplotype frequencies

Joins microhaplotype allele frequencies to translated amino acid calls,
aggregates frequencies per amino acid allele, and renormalises so
frequencies sum to one. Collapsed output averages evenly across
overlapping targets (taking `max(sample_total)`).

## Usage

``` r
slaf_from_mhaps_freqs(
  mhaps_slaf,
  loci_of_interest_per_microhaps,
  slaf_output = NULL,
  per_target_slaf_output = NULL
)
```

## Arguments

- mhaps_slaf:

  Path or data frame of microhaplotype SLAF. See *Inputs*.

- loci_of_interest_per_microhaps:

  Path or data frame of translated loci. See *Inputs*.

- slaf_output:

  Optional path for collapsed SLAF TSV.

- per_target_slaf_output:

  Optional path for per-target SLAF TSV.

## Value

A list with tibbles `slaf` (collapsed) and `per_target_slaf`.

## Details

### Inputs

- **`mhaps_slaf`**: Microhaplotype SLAF (`target_name`, `seq`, `freq`,
  `sample_total`), as a file path or data frame. See
  [`vignette("input-formats", package = "PGEcore")`](https://plasmogenepi.github.io/PGEcore/articles/input-formats.md).

- **`loci_of_interest_per_microhaps`**: Translated loci per
  microhaplotype (`target_name`, `gene_id`, `aa_position`, `seq`, `aa`),
  as a file path or data frame.

### Outputs

- **`slaf_output`** (optional): Collapsed SLAF TSV (`variant`, `freq`,
  `sample_total`).

- **`per_target_slaf_output`** (optional): Per-target SLAF TSV
  (`target_name`, `variant`, `freq`, `sample_total`).

### Running

    slaf_from_mhaps_freqs(
      mhaps_slaf = "mhaps_slaf.tsv",
      loci_of_interest_per_microhaps = "loci_of_interest_per_microhaps.tsv",
      slaf_output = "slaf.tsv"
    )

    Rscript exec/slaf_from_mhaps_freqs \
      --mhaps_slaf mhaps_slaf.tsv \
      --loci_of_interest_per_microhaps loci_of_interest_per_microhaps.tsv \
      --slaf_output slaf.tsv

## See also

[`vignette("input-formats", package = "PGEcore")`](https://plasmogenepi.github.io/PGEcore/articles/input-formats.md)
