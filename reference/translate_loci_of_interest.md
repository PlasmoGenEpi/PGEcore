# Translate loci of interest from microhaplotype sequences

Aligns each unique haplotype to its panel reference with an overlap
pairwise alignment, extracts the codon at each locus of interest, and
translates it.

## Usage

``` r
translate_loci_of_interest(
  allele_table,
  ref_bed,
  loci_of_interest,
  output_dir,
  select_target_names = NULL,
  select_specimen_names = NULL,
  overwrite_dir = FALSE,
  output_stop_codons = FALSE,
  collapse_calls_by_summing = FALSE
)
```

## Arguments

- allele_table:

  Path or data frame of allele table. See *Inputs*.

- ref_bed:

  Path or data frame of panel BED with `ref_seq`. See *Inputs*.

- loci_of_interest:

  Path or data frame of codon BED. See *Inputs*.

- output_dir:

  Directory to write results. Created if missing.

- select_target_names:

  Optional comma-separated names, path to a one-column TSV, or character
  vector of targets to keep.

- select_specimen_names:

  Optional comma-separated names, path to a one-column TSV, or character
  vector of specimens to keep.

- overwrite_dir:

  If `FALSE` (default), refuse to replace an existing `output_dir`.

- output_stop_codons:

  If `FALSE` (default), treat `*` as untranslatable (along with `X`).

- collapse_calls_by_summing:

  If `TRUE`, sum reads across overlapping targets; otherwise keep the
  target with the highest read count.

## Value

A named list with `loci_of_interest_for_target_for_microhap`,
`amino_acid_calls`, `collapsed_amino_acid_calls`,
`loci_covered_by_target_samples_info`, and `untranslatable`.

## Details

### Inputs

- **`allele_table`**: Allele table (`specimen_name`, `target_name`,
  `reads`, `seq`), as a file path or data frame. See
  [`vignette("input-formats", package = "PGEcore")`](https://plasmogenepi.github.io/PGEcore/articles/input-formats.md).

- **`ref_bed`**: Panel BED with `ref_seq` (`#chrom`, `start`, `end`,
  `target_name`, `length`, `strand`, `ref_seq`).

- **`loci_of_interest`**: Codon BED (`#chrom`, `start`, `end`, `name`,
  `length`, `strand`, `gene`, `gene_id`, `aa_position`); each locus must
  have `length == 3`.

### Outputs

- **`output_dir`**: Directory receiving
  `loci_of_interest_for_target_for_microhap.tsv.gz`,
  `amino_acid_calls.tsv.gz`, `collapsed_amino_acid_calls.tsv.gz`,
  `loci_covered_by_target_samples_info.tsv`, and optionally
  `allele_table_out_untranslatable.tsv`.

### Running

    translate_loci_of_interest(
      allele_table = "allele_table.tsv",
      ref_bed = "ref_bed_with_seq.tsv",
      loci_of_interest = "loci.bed",
      output_dir = "translate_out"
    )

    Rscript exec/translate_loci_of_interest \
      --allele_table allele_table.tsv \
      --ref_bed ref_bed_with_seq.tsv \
      --loci_of_interest loci.bed \
      --output_dir translate_out

Requires **Biostrings** and **pwalign** (Suggests).

## See also

[`vignette("input-formats", package = "PGEcore")`](https://plasmogenepi.github.io/PGEcore/articles/input-formats.md)
