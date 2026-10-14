# Split a SNP-Slice result into single-chain objects

Identical haplotypes are always collapsed first; see
[`snpslice_dedup_matrices()`](https://plasmogenepi.github.io/PGEcore/reference/snpslice_dedup_matrices.md).

## Usage

``` r
snpslice_split_chains(snpslice_res)
```

## Arguments

- snpslice_res:

  Result of
  [`snp.slicer::snp_slice()`](https://plasmogenepi.github.io/snp.slicer/reference/snp_slice.html).

## Value

List with `chains` (single-chain results in chain order) and `best`, the
index of the chain SNP-Slice reports (highest MAP log posterior).
