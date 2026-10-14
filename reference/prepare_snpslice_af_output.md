# Format SNP-Slice allele frequencies with variantstring names

`freq` is estimated from the choosen chain SNP-Slice reports (highest
MAP log posterior by default) where as `averaged_freq` is similar to how
Ju et al. (2024) did COI with taking the average across the chains.
`averaged_freq` is the mean over all chains's frequency for the
haplotype, with a haplotype a chain never produced counting as zero
there. A row is kept when either column is non-zero This increases
sensitivity but also decrease precision

## Usage

``` r
prepare_snpslice_af_output(chains, best, loci_groups, estimator)
```

## Arguments

- chains, best:

  Output of
  [`snpslice_split_chains()`](https://plasmogenepi.github.io/PGEcore/reference/snpslice_split_chains.md).
