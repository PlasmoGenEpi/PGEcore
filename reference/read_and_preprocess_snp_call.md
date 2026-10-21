# Read SNP call data for THEREALMcCOIL

Validates the required columns and keeps only biallelic loci with an
observed reference base (see
[`filter_biallelic()`](https://plasmogenepi.github.io/PGEcore/reference/filter_biallelic.md)).

## Usage

``` r
read_and_preprocess_snp_call(snp_calls)
```

## Arguments

- snp_calls:

  Path to an independent, collapsed SNP call TSV.

## Value

A data frame with `specimen_name`, `snp_name`, `ref_base`, `seq_base`,
`reads`.
