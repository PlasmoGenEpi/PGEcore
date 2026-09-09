# Estimate IBD-based relatedness with Dcifer

Estimates pairwise relatedness from an allele table (optionally with
COI, allele frequencies, and population metadata). Requires **dcifer**
and parallel helpers (**doParallel**, **parallelly**, **foreach**,
**iterators**) (Suggests).

## Usage

``` r
dcifer_ibd_wrapper(
  allele_table,
  relatedness_output,
  coi_table = NULL,
  allele_freq_table = NULL,
  specimen_name_col = "specimen_name",
  target_name_col = "target_name",
  target_value_col = "seq",
  specimen_metadata = NULL,
  pop_name_col = NULL,
  rnull = 0,
  alpha = 0.05,
  use_estm = FALSE,
  threads = 1L,
  verbose = FALSE
)
```

## Arguments

- allele_table:

  Path to allele table TSV. See *Inputs*.

- relatedness_output:

  Path for relatedness TSV. See *Outputs*.

- coi_table:

  Optional path to COI table TSV. See *Inputs*.

- allele_freq_table:

  Optional path to allele-frequency TSV. See *Inputs*.

- specimen_name_col, target_name_col, target_value_col:

  Column names in the allele / COI / frequency tables.

- specimen_metadata:

  Optional metadata TSV. See *Inputs*.

- pop_name_col:

  Optional population column in metadata.

- rnull:

  Relatedness null for hypothesis testing.

- alpha:

  Significance level.

- use_estm:

  If `TRUE`, use
  [`dcifer::ibdEstM()`](https://eppicenter.github.io/dcifer/reference/ibdEstM.html)
  instead of
  [`dcifer::ibdPair()`](https://eppicenter.github.io/dcifer/reference/ibdPair.html).

- threads:

  Number of parallel workers.

- verbose:

  Print parallel worker output.

## Value

The relatedness tibble (also written to `relatedness_output`).

## Details

### Inputs

- **`allele_table`**: Allele table TSV. See
  [`vignette("input-formats", package = "PGEcore")`](https://plasmogenepi.github.io/PGEcore/articles/input-formats.md).

- **`coi_table`**: Optional COI table TSV. If omitted, COI is inferred.

- **`allele_freq_table`**: Optional SLAF TSV (`target_name`, `seq`,
  `freq`). If omitted, frequencies are estimated from the allele table.

- **`specimen_metadata`**: Optional metadata TSV for population-specific
  runs (with `pop_name_col`).

### Outputs

- **`relatedness_output`**: Relatedness TSV with `specimen_name_a`,
  `specimen_name_b`, `btwn_host_rel`, and optionally `p_value`,
  `CI_lower`, `CI_upper`, and/or `strain_pair`.

### Running

    dcifer_ibd_wrapper(
      allele_table = "allele_table.tsv",
      relatedness_output = "relatedness.tsv"
    )

    Rscript exec/dcifer_ibd_wrapper \
      --allele_table allele_table.tsv \
      --relatedness_output relatedness.tsv

Requires **dcifer** and parallel Suggests packages.

## See also

[`vignette("input-formats", package = "PGEcore")`](https://plasmogenepi.github.io/PGEcore/articles/input-formats.md),
[`dcifer_slaf_wrapper()`](https://plasmogenepi.github.io/PGEcore/reference/dcifer_slaf_wrapper.md)
