# Run IDM/OM MLE independently at each locus

Run IDM/OM MLE independently at each locus

## Usage

``` r
run_idm_mle_across_loci(
  df,
  model = "IDM",
  lambda_initial = 1,
  eps_initial = 0.1
)
```

## Arguments

- df:

  Formatted input from
  [`prepare_input_4_allele_table()`](https://plasmogenepi.github.io/PGEcore/reference/prepare_input_4_allele_table.md)
  or
  [`prepare_input_4_aa_calls()`](https://plasmogenepi.github.io/PGEcore/reference/prepare_input_4_aa_calls.md).

- model:

  `"IDM"`, `"OM"`, `"IDM_OM"`, `"OM_CC"`, or `"IDM_OM_CC"`.

- lambda_initial:

  Initial lambda for the numerical solver.

- eps_initial:

  Initial epsilon for the numerical solver.

## Value

Tibble with `variant` and `freq`.
