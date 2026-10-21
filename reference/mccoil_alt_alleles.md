# Reference and alternate allele of each biallelic SNP

Both McCOIL models estimate the frequency of the non-reference allele.

## Usage

``` r
mccoil_alt_alleles(df)
```

## Arguments

- df:

  Output of
  [`read_and_preprocess_snp_call()`](https://plasmogenepi.github.io/PGEcore/reference/read_and_preprocess_snp_call.md).

## Value

A tibble with `snp_name`, `ref_base` and `alt_base`.
