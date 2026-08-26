# Compute MCMC convergence diagnostics across McCOIL chains

Assembles a posterior draws array from the post-burn-in portion of the
per-chain traces and summarises it with
[`summarize_convergence_draws()`](https://plasmogenepi.github.io/PGEcore/reference/summarize_convergence_draws.md).

## Usage

``` r
prepare_mccoil_convergence_output(traces, summary_path, totalrun, burnin)
```

## Arguments

- traces:

  Per-chain trace data frames from
  [`run_mccoil_chains()`](https://plasmogenepi.github.io/PGEcore/reference/run_mccoil_chains.md).

- summary_path:

  Path to the `*_summary.txt` written by chain 1, used for the parameter
  names and their order in the traces.

- totalrun:

  Total MCMC iterations.

- burnin:

  Burn-in iterations.

## Value

A data frame of convergence diagnostics, one row per parameter.
