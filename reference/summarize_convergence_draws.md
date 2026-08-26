# Summarize posterior draws with convergence diagnostics

Keeps the subset of
[`posterior::summarise_draws()`](https://mc-stan.org/posterior/reference/draws_summary.html)
columns reported by the MCMC wrappers, so convergence tables share one
schema across tools.

## Usage

``` r
summarize_convergence_draws(draws)
```

## Arguments

- draws:

  A `draws` object (for example from
  [`posterior::as_draws_array()`](https://mc-stan.org/posterior/reference/draws_array.html)).

## Value

A data frame with one row per parameter and the columns `variable`,
`mean`, `median`, `sd`, `q5`, `q95`, `rhat`, `ess_bulk`, and `ess_tail`.
