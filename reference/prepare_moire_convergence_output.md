# Compute MCMC convergence diagnostics across all chains

Assembles a posterior draws array (iteration x chain x variable)
covering every estimated parameter (per-sample COI,
false-positive/false-negative error rates, within-host relatedness;
per-locus/allele frequencies; and the population mean COI) and
summarizes it with
[`summarize_convergence_draws()`](https://plasmogenepi.github.io/PGEcore/reference/summarize_convergence_draws.md).

## Usage

``` r
prepare_moire_convergence_output(mcmc_results)
```

## Arguments

- mcmc_results:

  The list returned by
  [`run_moire()`](https://plasmogenepi.github.io/PGEcore/reference/run_moire.md).

## Value

A data frame of convergence diagnostics, one row per parameter.

## Details

Chains are not guaranteed to be the same length: MOIRE's `max_runtime`
stops each chain independently once its own wall clock expires, so a
truncated run yields ragged chains. All chains are truncated to the
shortest one (with a warning) so the draws array is rectangular and
iteration `i` refers to the same sweep in every chain.
