# Does the installed MOIRE accept a seed?

Feature detection rather than a version comparison: whether
[`moire::run_mcmc()`](https://EPPIcenter.github.io/moire/reference/run_mcmc.html)
takes a `seed` argument is exactly the question, and the formals answer
it directly. A version string would need updating every time the feature
moves between branches or releases.

## Usage

``` r
moire_supports_seed()
```

## Value

`TRUE` when
[`moire::run_mcmc()`](https://EPPIcenter.github.io/moire/reference/run_mcmc.html)
has a `seed` argument.
