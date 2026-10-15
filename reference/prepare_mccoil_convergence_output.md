# Compute MCMC convergence diagnostics across McCOIL chains

Compute MCMC convergence diagnostics across McCOIL chains

## Usage

``` r
prepare_mccoil_convergence_output(draws)
```

## Arguments

- draws:

  Output of
  [`mccoil_draws()`](https://plasmogenepi.github.io/PGEcore/reference/mccoil_draws.md).

## Value

A data frame of convergence diagnostics, one row per parameter.
