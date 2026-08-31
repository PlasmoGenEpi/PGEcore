# Stop early when MOIRE recorded no draws for some chain

MOIRE's `max_runtime` stops each chain the moment its own wall clock
expires, which for a short enough cap happens partway through burn-in –
leaving that chain with no recorded draws at all. Nothing downstream
copes with that: moire's own summarizers fail deep inside
[`quantile()`](https://rdrr.io/r/stats/quantile.html) with "'x' must be
atomic", and
[`extract_moire_chain_draws()`](https://plasmogenepi.github.io/PGEcore/reference/extract_moire_chain_draws.md)
fails reading the allele count off a first draw that does not exist.
Both are opaque, so the condition is caught here instead.

## Usage

``` r
assert_moire_chains_have_draws(mcmc_results)
```

## Arguments

- mcmc_results:

  The list returned by
  [`run_moire()`](https://plasmogenepi.github.io/PGEcore/reference/run_moire.md).

## Value

`invisible(NULL)`; called for its side effect of erroring.
