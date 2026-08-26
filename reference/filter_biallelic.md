# Restrict SNP calls to biallelic loci for the proportional model

The proportional model represents each locus with exactly two alleles
(`a1`/`a2`), so loci with any other number of alleles cannot be encoded
correctly. Loci with exactly two distinct `seq_base` values are kept;
monomorphic and multiallelic loci are dropped with a warning.

## Usage

``` r
filter_biallelic(df)
```

## Arguments

- df:

  Output of
  [`read_and_preprocess_snp_call()`](https://plasmogenepi.github.io/PGEcore/reference/read_and_preprocess_snp_call.md).

## Value

`df` filtered to biallelic loci only.
