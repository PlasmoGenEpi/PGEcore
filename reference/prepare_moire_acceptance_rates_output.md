# Compute parallel tempering swap acceptance rates

For each independent chain and each temperature rung, computes the swap
(exchange) acceptance rate with the adjacent hotter rung, following
MOIRe's own convention (see
[`moire::plot_chain_swaps()`](https://EPPIcenter.github.io/moire/reference/plot_chain_swaps.html)):
`swap_acceptances / (samples_per_chain / 2)`.

## Usage

``` r
prepare_moire_acceptance_rates_output(mcmc_results)
```

## Arguments

- mcmc_results:

  The list returned by
  [`run_moire()`](https://plasmogenepi.github.io/PGEcore/reference/run_moire.md).

## Value

A tibble with one row per chain-rung combination and the columns
`chain`, `rung`, `temperature`, and `swap_acceptance_rate`.

## Details

Swap acceptances are recorded per adjacent rung pair, so the rate for
rung `k` describes swaps between rung `k` and rung `k + 1`; the final
rung has no partner above it and its rate is `NA`. `temperature` is
MOIRe's `temp_gradient` value for the rung (the power-posterior exponent
in `[0, 1]`), read per chain as it may be adapted.
