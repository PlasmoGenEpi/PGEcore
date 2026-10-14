# Apply [`snpslice_dedup_matrices()`](https://plasmogenepi.github.io/PGEcore/reference/snpslice_dedup_matrices.md) to every estimate a chain carries

Apply
[`snpslice_dedup_matrices()`](https://plasmogenepi.github.io/PGEcore/reference/snpslice_dedup_matrices.md)
to every estimate a chain carries

## Usage

``` r
snpslice_dedup_chain(chain)
```

## Arguments

- chain:

  A single-chain `snp_slice_results` object from
  [`snp.slicer::get_chain()`](https://plasmogenepi.github.io/snp.slicer/reference/get_chain.html).

## Value

The chain with its MAP, final-sample and (if present) stored MCMC sample
matrices de-duplicated.
