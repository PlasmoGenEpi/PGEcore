# Format pooled McCOIL draws into PGE COI and SLAF tables

Estimates are posterior medians over the post-burn-in draws of all
chains. Each SLAF row is the non-reference allele of one SNP, named
`snp_name:seq_base` in `variant`.

## Usage

``` r
format_mccoil_output(draws, alleles)
```

## Arguments

- draws:

  Output of
  [`mccoil_draws()`](https://plasmogenepi.github.io/PGEcore/reference/mccoil_draws.md).

- alleles:

  Output of
  [`mccoil_alt_alleles()`](https://plasmogenepi.github.io/PGEcore/reference/mccoil_alt_alleles.md).

## Value

A list with `slaf` (`variant`, `snp_name`, `ref_base`, `seq_base`,
`freq`) and `coi` tibbles.
