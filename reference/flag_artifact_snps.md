# Flag SNPs whose minor allele recurs in otherwise-monoclonal specimens

A real minor strain shows minor alleles at many SNPs of a specimen, so a
SNP whose minor allele keeps appearing in specimens with at most one
minor allele elsewhere is an artifact (e.g. PCR stutter in a tandem
repeat).

## Usage

``` r
flag_artifact_snps(
  snp_table,
  min_minor_rate = 0.02,
  min_samples = 3,
  min_depth = 30
)
```

## Arguments

- snp_table:

  SNP call tibble.

- min_minor_rate:

  Minimum share of eligible specimens showing the SNP's minor allele.

- min_samples:

  Minimum number of eligible specimens showing it.

- min_depth:

  Minimum reads at a SNP for a specimen to be considered.

## Value

A tibble with one row per SNP: `snp_name`, `eligible_specimens`
(specimens with at most one minor allele at other SNPs),
`total_minor_allele_detected` (those with a minor allele at this SNP),
`minor_rate` and `artifact`.
