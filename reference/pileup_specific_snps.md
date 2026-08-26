# Pile up specific SNPs covered by microhaplotype sequences

Aligns each unique haplotype to its panel reference with an overlap
pairwise alignment, then extracts bases at SNP-of-interest coordinates.

## Usage

``` r
pileup_specific_snps(
  allele_table,
  ref_bed,
  snps_of_interest,
  output_dir,
  select_target_names = NULL,
  select_specimen_names = NULL,
  overwrite_dir = FALSE,
  collapse_calls_by_summing = FALSE
)
```

## Arguments

- allele_table:

  Path or data frame of allele table. See *Inputs*.

- ref_bed:

  Path or data frame of panel BED with `ref_seq`. See *Inputs*.

- snps_of_interest:

  Path or data frame of SNP BED. See *Inputs*.

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

- collapse_calls_by_summing:

  If `TRUE`, sum reads across overlapping targets; otherwise keep the
  target with the highest read count.

## Value

A named list with `snp_calls`, `collapsed_snp_calls`,
`snps_covered_by_target_samples_info`, and `uncallable`.

## Details

### Inputs

- **`allele_table`**: Allele table (`specimen_name`, `target_name`,
  `reads`, `seq`), as a file path or data frame. See
  [`vignette("input-formats", package = "PGEcore")`](https://plasmogenepi.github.io/PGEcore/articles/input-formats.md).

- **`ref_bed`**: Panel BED with `ref_seq` (`#chrom`, `start`, `end`,
  `target_name`, `length`, `strand`, `ref_seq`).

- **`snps_of_interest`**: SNP BED (`#chrom`, `start`, `end`, `name`,
  `length`, `strand`); each SNP must span one base (`end - start == 1`).

### Outputs

- **`output_dir`**: Directory receiving `snp_calls.tsv.gz`,
  `collapsed_snp_calls.tsv.gz`,
  `snps_covered_by_target_samples_info.tsv`, and optionally
  `allele_table_out_uncallable.tsv`.

### Running

    pileup_specific_snps(
      allele_table = "allele_table.tsv",
      ref_bed = "ref_bed_with_seq.tsv",
      snps_of_interest = "snps.bed",
      output_dir = "pileup_out"
    )

    Rscript exec/pileup_specific_snps \
      --allele_table allele_table.tsv \
      --ref_bed ref_bed_with_seq.tsv \
      --snps_of_interest snps.bed \
      --output_dir pileup_out

Requires **Biostrings** and **pwalign** (Suggests).

## See also

[`vignette("input-formats", package = "PGEcore")`](https://plasmogenepi.github.io/PGEcore/articles/input-formats.md)
