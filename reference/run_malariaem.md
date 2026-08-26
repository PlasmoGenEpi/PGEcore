# Run malaria.em and write frequency and phase summaries

Runs malaria.em haplotype EM on an in-memory allele matrix and writes
genotype-frequency and phasing summaries. Requires **malaria.em**
(Suggests). For reading allele tables from disk, use
[`malariaem_wrapper()`](https://plasmogenepi.github.io/PGEcore/reference/malariaem_wrapper.md).

## Usage

``` r
run_malariaem(
  matrix,
  test_size = "min",
  max_size = "8",
  subset_targets = FALSE,
  target_groups = NULL,
  freq_output = NULL,
  phase_output = NULL
)
```

## Arguments

- matrix:

  Allele matrix (specimens as rows, loci as columns). See *Inputs*.

- test_size:

  Maximum COI to test, or `"min"` for the minimum allowed.

- max_size:

  Error if inferred COI range exceeds this cutoff.

- subset_targets:

  If `TRUE`, run separately for each `group_id`.

- target_groups:

  Data frame with `group_id` and `target_name`. See *Inputs*.

- freq_output:

  Path for genotype-frequency TSV. See *Outputs*.

- phase_output:

  Path for phasing TSV. See *Outputs*.

## Value

A list of malaria.em results (or a named list per group).

## Details

### Inputs

- **`matrix`**: Allele matrix (specimens as rows, loci as columns;
  alleles space-separated within cells).

- **`target_groups`**: Optional data frame with `group_id` and
  `target_name` when `subset_targets = TRUE`.

### Outputs

- **`freq_output`**: Genotype-frequency TSV (`gt_id`, `target_name`,
  `seq`, `freq`, `freq_se`; plus `group_id` when subsetting).

- **`phase_output`**: Phasing TSV (`specimen_name`, `target_name`,
  `seq`, `gt_id`, `posterior_est`, `phase_id`; plus `group_id` when
  subsetting).

### Running

    run_malariaem(
      matrix = allele_matrix,
      freq_output = "gt_freq_summary_all.tsv",
      phase_output = "gt_phase_summary_all.tsv"
    )

File and CLI users should call
[`malariaem_wrapper()`](https://plasmogenepi.github.io/PGEcore/reference/malariaem_wrapper.md)
/ `Rscript exec/malariaem_wrapper ...`.

Requires **malaria.em** (Suggests).

## See also

[`malariaem_wrapper()`](https://plasmogenepi.github.io/PGEcore/reference/malariaem_wrapper.md),
[`vignette("input-formats", package = "PGEcore")`](https://plasmogenepi.github.io/PGEcore/articles/input-formats.md)
