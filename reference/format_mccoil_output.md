# Format pooled McCOIL draws into PGE COI and SLAF tables

Estimates are posterior medians over the post-burn-in draws of all
chains.

## Usage

``` r
format_mccoil_output(draws)
```

## Arguments

- draws:

  Output of
  [`mccoil_draws()`](https://plasmogenepi.github.io/PGEcore/reference/mccoil_draws.md).

## Value

A list with `slaf` and `coi` tibbles.
