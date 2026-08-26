# Run malaria.em from allele-table and output paths

Reads an allele table, builds the EM matrix, and calls
[`run_malariaem()`](https://plasmogenepi.github.io/PGEcore/reference/run_malariaem.md).
Requires **malaria.em** (Suggests).

## Usage

``` r
malariaem_wrapper(
  allele_table,
  subset_targets = FALSE,
  target_groups = NULL,
  max_size = "8",
  test_size = "min",
  freq_output = "gt_freq_summary_all.tsv",
  phase_output = "gt_phase_summary_all.tsv",
  seed = 1L
)
```

## Arguments

- allele_table:

  Path to allele table TSV. See *Inputs*.

- subset_targets:

  Logical; subset by `target_groups`.

- target_groups:

  Optional path to groups TSV. See *Inputs*.

- max_size:

  COI cutoff (legacy CLI default `"8"`).

- test_size:

  COI size to test, or `"min"`.

- freq_output:

  Frequency summary path. See *Outputs*.

- phase_output:

  Phase summary path. See *Outputs*.

- seed:

  Random seed.

## Value

The object returned by
[`run_malariaem()`](https://plasmogenepi.github.io/PGEcore/reference/run_malariaem.md).

## Details

### Inputs

- **`allele_table`**: Allele table TSV (`specimen_name`, `target_name`,
  `seq`). See
  [`vignette("input-formats", package = "PGEcore")`](https://plasmogenepi.github.io/PGEcore/articles/input-formats.md).

- **`target_groups`**: Optional groups TSV (`group_id`, `target_name`)
  when `subset_targets = TRUE`.

### Outputs

- **`freq_output`**: Genotype-frequency TSV (`gt_id`, `target_name`,
  `seq`, `freq`, `freq_se`).

- **`phase_output`**: Phasing TSV (`specimen_name`, `target_name`,
  `seq`, `gt_id`, `posterior_est`, `phase_id`).

### Running

    malariaem_wrapper(
      allele_table = "allele_table.tsv",
      freq_output = "gt_freq_summary_all.tsv",
      phase_output = "gt_phase_summary_all.tsv"
    )

    Rscript exec/malariaem_wrapper \
      --allele_table allele_table.tsv \
      --freq_output gt_freq_summary_all.tsv \
      --phase_output gt_phase_summary_all.tsv

Requires **malaria.em** (Suggests).

## See also

[`run_malariaem()`](https://plasmogenepi.github.io/PGEcore/reference/run_malariaem.md),
[`vignette("input-formats", package = "PGEcore")`](https://plasmogenepi.github.io/PGEcore/articles/input-formats.md)
