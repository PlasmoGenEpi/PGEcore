# Format SNP-Slice COI estimates

Three per-specimen columns:

- `coi`: number of strains assigned in the reported chain (highest MAP
  log posterior), the row sum of its allocation matrix.

- `coi_cons_weighted`:
  [`snpslice_consensus_coi()`](https://plasmogenepi.github.io/PGEcore/reference/snpslice_consensus_coi.md)
  across chains.

- `coi_chain_mean`: the mean over all chains of that chain's `coi`,
  which is how Ju et al. (2024) report COI (the average of the
  final-sample estimate over independently initialised chains).

## Usage

``` r
prepare_snpslice_coi_output(chains, best, specimen_name_col, estimator)
```

## Arguments

- chains, best:

  Output of
  [`snpslice_split_chains()`](https://plasmogenepi.github.io/PGEcore/reference/snpslice_split_chains.md).
