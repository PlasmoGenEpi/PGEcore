# Restrict SNP calls to biallelic loci with an observed reference base

Both McCOIL models represent a locus as a reference and one alternate
allele. Loci are kept when exactly two distinct `seq_base` values are
observed and one of them is `ref_base`; monomorphic and multiallelic
loci, and loci whose reference base is never observed, are dropped with
a warning.

## Usage

``` r
filter_biallelic(df)
```

## Arguments

- df:

  SNP calls with `snp_name`, `ref_base`, `seq_base` and `reads`.

## Value

`df` filtered to those loci.
