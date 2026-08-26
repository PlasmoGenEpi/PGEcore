# Filter SNPs to highest-diversity loci spaced by a minimum distance

Ranks SNPs by expected heterozygosity, then greedily keeps loci that are
at least `mindist_between_snps` apart on the same chromosome.

## Usage

``` r
filter_to_highest_diversity_independent_snp_call(
  snp_calls,
  snp_calls_output = NULL,
  mindist_between_snps = 10000,
  select_target_names = NULL,
  select_specimen_names = NULL,
  overwrite = FALSE,
  only_biallelic = FALSE,
  only_informative = FALSE
)
```

## Arguments

- snp_calls:

  Path to a SNP calls TSV, or a data frame with the same columns. See
  *Inputs*.

- snp_calls_output:

  Optional output TSV path.

- mindist_between_snps:

  Minimum distance between kept SNPs (default `10000`).

- select_target_names:

  Optional comma-separated names, path to a one-column TSV, or character
  vector of targets to keep.

- select_specimen_names:

  Optional comma-separated names, path to a one-column TSV, or character
  vector of specimens to keep.

- overwrite:

  If `FALSE` (default), refuse to overwrite `snp_calls_output`.

- only_biallelic:

  If `TRUE`, restrict to rows where `is_biallelic` is `TRUE`.

- only_informative:

  If `TRUE`, drop SNPs with expected heterozygosity of zero.

## Value

Filtered SNP table including an `he` column.

## Details

### Inputs

- **`snp_calls`**: SNP calls (`specimen_name`, `target_name`, `chrom`,
  `pos`, `snp_name`, `ref_base`, `seq_base`, `reads`, `is_biallelic`),
  as a file path or data frame. See
  [`vignette("input-formats", package = "PGEcore")`](https://plasmogenepi.github.io/PGEcore/articles/input-formats.md).

### Outputs

- **`snp_calls_output`** (optional): Filtered SNP calls TSV including an
  `he` column. If `NULL`, results are returned without writing a file.

### Running

    filter_to_highest_diversity_independent_snp_call(
      snp_calls = "snp_calls.tsv",
      snp_calls_output = "filtered_snp_calls.tsv"
    )

    Rscript exec/filter_to_highest_diversity_independent_snp_call \
      --snp_calls snp_calls.tsv \
      --snp_calls_output filtered_snp_calls.tsv

## See also

[`vignette("input-formats", package = "PGEcore")`](https://plasmogenepi.github.io/PGEcore/articles/input-formats.md)
