# Write IDM SLAF output

Write IDM SLAF output

## Usage

``` r
write_idm_output(res, slaf_output, allele_table = FALSE)
```

## Arguments

- res:

  Result table from
  [`run_idm_mle_across_loci()`](https://plasmogenepi.github.io/PGEcore/reference/run_idm_mle_across_loci.md).

- slaf_output:

  Output TSV path.

- allele_table:

  If `TRUE`, split `variant` into `target_name` and `seq`.
