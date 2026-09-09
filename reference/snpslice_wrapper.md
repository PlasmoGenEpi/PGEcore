# Estimate multilocus allele frequency and COI with SNP-Slice

Estimates multilocus allele frequencies and per-specimen COI. Requires
**snp.slicer** (with multi-chain sampling and `estimate`) and
**variantstring** 1.x (Suggests).

## Usage

``` r
snpslice_wrapper(
  allele_table,
  loci_groups,
  mlaf_output,
  coi_output,
  convergence_output = "convergence_diag.tsv",
  specimen_name_col = "specimen_name",
  target_name_col = "aa_locus",
  target_value_col = "aa",
  target_count_col = "reads",
  loci_limit = NULL,
  model = "negative_binomial",
  n_sample = 10000L,
  n_burnin = NULL,
  alpha = 2.6,
  rho = NULL,
  threshold = 0.001,
  gap = NULL,
  estimator = "final_sample",
  n_chains = 3L,
  threads = 1L,
  verbose = FALSE,
  seed = 1L
)
```

## Arguments

- allele_table:

  Path to allele / AA-calls TSV with counts. See *Inputs*.

- loci_groups:

  Path to loci-groups TSV. See *Inputs*.

- mlaf_output:

  Path for multilocus allele-frequency TSV. See *Outputs*.

- coi_output:

  Path for COI TSV. See *Outputs*.

- convergence_output:

  Path for per-restart optimization-diagnostics TSV. See *Outputs*.

- specimen_name_col, target_name_col, target_value_col,
  target_count_col:

  Column names in `allele_table`.

- loci_limit:

  Optional cap on the number of loci.

- model:

  Observation model for SNP-Slice.

- n_sample:

  Post-burn-in MCMC iterations retained per chain.

- n_burnin:

  Burn-in iterations per chain. If `NULL`, SNP-Slice uses
  `floor(n_sample / 2)`.

- alpha, threshold, gap:

  SNP-Slice MCMC settings.

- rho:

  Dictionary sparsity parameter. If `NULL`, it is left to
  [`snp.slicer::snp_slice()`](https://plasmogenepi.github.io/snp.slicer/reference/snp_slice.html),
  which defaults to 0.5 for the categorical model and the global minor
  allele frequency for count models.

- estimator:

  Estimator for COI and allele frequencies: `"final_sample"` (default,
  matching the SNP-Slice paper), `"map"`, or `"posterior"` (the
  posterior mean, which also yields uncertainty columns).

- n_chains:

  Number of independent MCMC chains. More than one is required for the
  Gelman-Rubin R-hat diagnostic.

- threads:

  Cores used to run chains simultaneously (capped at `n_chains`).

- verbose:

  Verbose SNP-Slice output.

- seed:

  Random seed.

## Value

Invisibly, a list with `mlaf`, `coi`, and `convergence` tibbles.

## Details

### Inputs

- **`allele_table`**: Allele / AA-style table with counts. Default
  column names map AA-call fields (`aa_locus`, `aa`, `reads`). See
  [`vignette("input-formats", package = "PGEcore")`](https://plasmogenepi.github.io/PGEcore/articles/input-formats.md).

- **`loci_groups`**: Loci-groups TSV (`group_id` plus a locus column
  matching `target_name_col`).

### Outputs

- **`mlaf_output`**: Multilocus allele frequencies (`group_id`,
  `variant`, `freq`, …).

- **`coi_output`**: COI estimates (`specimen_name`, `coi`,
  `coi_cons_weighted`; uncertainty columns when
  `estimator = "posterior"`). `coi` counts every strain assigned to a
  host in the best restart. `coi_cons_weighted` pools haplotype
  membership across all restarts and weights each haplotype by its
  consensus support, which counters the dictionary over-parameterisation
  that inflates `coi`.

- **`convergence_output`**: Per-restart optimization diagnostics
  (`chain_id`, `seed`, `map_logpost`, `is_best`, `map_iteration`,
  `final_iteration`, `plateau_frac`, `map_kstar`, `map_ktrunc`,
  `coi_mean`, `coi_ccc_to_best`). SNP-Slice reports the restart with the
  highest MAP log posterior rather than pooling chains, so between-chain
  R-hat and ESS do not describe its output and are not emitted.

### Running

    snpslice_wrapper(
      allele_table = "aa_calls.tsv",
      loci_groups = "loci_groups.tsv",
      mlaf_output = "mlaf.tsv",
      coi_output = "coi.tsv"
    )

    Rscript exec/snpslice_wrapper \
      --allele_table aa_calls.tsv \
      --loci_groups loci_groups.tsv \
      --mlaf_output mlaf.tsv \
      --coi_output coi.tsv

Requires **snp.slicer** and **variantstring** (Suggests).

## See also

[`vignette("input-formats", package = "PGEcore")`](https://plasmogenepi.github.io/PGEcore/articles/input-formats.md)
