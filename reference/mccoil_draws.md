# Assemble post-burn-in McCOIL draws across chains

Assemble post-burn-in McCOIL draws across chains

## Usage

``` r
mccoil_draws(traces, summary_path, totalrun, burnin)
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

A list with `draws` (iteration x chain x variable array) and the
`coi_names` and `freq_names` it covers.
