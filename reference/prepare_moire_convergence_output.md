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
