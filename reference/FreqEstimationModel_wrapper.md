# Estimate multilocus allele frequencies with FreqEstimationModel

Estimates multilocus haplotype frequencies from amino-acid calls and
COI. Requires **FreqEstimationModel**, **variantstring**, **posterior**,
and parallel helpers (**foreach**, **doMC**, **plyr**, **coda**,
**abind**) (Suggests).

## Usage

``` r
FreqEstimationModel_wrapper(
  aa_calls,
  coi,
  loci_groups,
  mlaf_output,
  threads = 1L,
  seed = 1L,
  n_chains = 3L,
  no_traces_preburnin = 10000L,
  thinning_interval = 1L,
  convergence_output = "convergence_diag.tsv"
)
```

## Arguments

- aa_calls:

  Path to amino-acid call TSV. See *Inputs*.

- coi:

  Path to COI table TSV, or a numeric average COI. See *Inputs*.

- loci_groups:

  Path to loci-groups TSV. See *Inputs*.

- mlaf_output:

  Output TSV path. See *Outputs*.

- threads:

  Number of threads.

- seed:

  Random seed.

- n_chains:

  Number of MCMC chains to run per group. At least two are needed to
  compute the Gelman-Rubin R-hat convergence diagnostic.

- no_traces_preburnin:

  Number of MCMC traces retained per chain before burn-in. Lower values
  cut memory roughly proportionally; see
  [`run_FreqEstimationModel()`](https://plasmogenepi.github.io/PGEcore/reference/run_FreqEstimationModel.md).

- thinning_interval:

  Metropolis-Hastings updates per retained trace.
  `no_traces_preburnin * thinning_interval` is the total sampler effort,
  so e.g. 2500 x 4 samples as far as the default 10000 x 1 while storing
  a quarter of the draws.

- convergence_output:

  Output TSV path for per-group MCMC convergence diagnostics. See
  *Outputs*.

## Value

The formatted output data frame (also written to `mlaf_output`).

## Details

### Inputs

- **`aa_calls`**: Amino-acid calls TSV. See
  [`vignette("input-formats", package = "PGEcore")`](https://plasmogenepi.github.io/PGEcore/articles/input-formats.md).

- **`coi`**: Path to a COI table TSV, **or** a numeric average COI.

- **`loci_groups`**: Loci-groups TSV (`group_id`, `gene_id`,
  `aa_position`).

### Outputs

- **`mlaf_output`**: Multilocus allele frequencies (`variant`, `freq`,
  `median_freq`, `CI_2.5`, `CI_97.5`, `sample_total`, `group_id`, …).

- **`convergence_output`**: Per-group MCMC diagnostics (`group_id`,
  `variable`, `mean`, `median`, `sd`, `q5`, `q95`, `rhat`, `ess_bulk`,
  `ess_tail`).

### Running

    FreqEstimationModel_wrapper(
      aa_calls = "aa_calls.tsv",
      coi = "coi_table.tsv",
      loci_groups = "loci_groups.tsv",
      mlaf_output = "mlaf.tsv"
    )

    Rscript exec/FreqEstimationModel_wrapper \
      --aa_calls aa_calls.tsv \
      --coi coi_table.tsv \
      --loci_groups loci_groups.tsv \
      --mlaf_output mlaf.tsv

Requires **FreqEstimationModel**, **variantstring**, **posterior**, and
parallel Suggests packages.

## See also

[`vignette("input-formats", package = "PGEcore")`](https://plasmogenepi.github.io/PGEcore/articles/input-formats.md)
