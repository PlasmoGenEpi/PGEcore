# Solve one locus with the vendored MLE

Solve one locus with the vendored MLE

## Usage

``` r
idm_locus_mle(tbl, model, lambda_initial, eps_initial, continuity = 0)
```

## Arguments

- tbl:

  Two-column table (`specimen_name`, `variants`) for a single locus.

- model:

  `"IDM"` or `"OM"` — the value handed to the vendored `MLE()`.

- lambda_initial:

  Initial lambda for the numerical solver.

- eps_initial:

  Initial epsilon for the numerical solver.

- continuity:

  Continuity correction applied to lineage prevalence; `0` is the
  published estimator.

## Value

List with `variants` and `freq`, one entry per variant at the locus.
