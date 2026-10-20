# Seed arguments to forward to `mcmc_sampling_parallel()`

`mcmc_sampling_parallel()` builds its per-chain RNG streams with
`RNGseq()`, so the seed has to be passed to it directly. When `seed` is
`NULL` a random seed is drawn; the wrapper records it in the convergence
summary.

## Usage

``` r
fem_seed_args(seed)
```

## Arguments

- seed:

  Integer seed, or `NULL` to draw one at random.

## Value

A list to splice into the `mcmc_sampling_parallel()` call: empty when
the installed build has no `seed` argument, otherwise
`list(seed = <int>)`.
