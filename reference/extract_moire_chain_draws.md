# Assemble a named list of every estimated parameter's draws for one chain

Assemble a named list of every estimated parameter's draws for one chain

## Usage

``` r
extract_moire_chain_draws(chain, sample_ids, loci)
```

## Arguments

- chain:

  One element of `mcmc_results$chains`.

- sample_ids:

  Character vector of specimen IDs; order matches the per-sample draw
  lists (`chain$coi`, `chain$eps_pos`, etc.).

- loci:

  Character vector of locus names; order matches `chain$allele_freqs`.

## Value

A named list of numeric draw vectors, each of length
`samples_per_chain`, named for the parameter it belongs to.
