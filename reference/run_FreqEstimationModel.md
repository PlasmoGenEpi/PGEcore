# Run FreqEstimationModel MCMC for one group

Run FreqEstimationModel MCMC for one group

## Usage

``` r
run_FreqEstimationModel(
  sample_matrix_list,
  COI,
  threads,
  seed,
  n_chains = 3L,
  no_traces_preburnin = 10000L,
  thinning_interval = 1L,
  moi_max = 8L
)
```

## Arguments

- sample_matrix_list:

  Output of
  [`create_FEM_input()`](https://plasmogenepi.github.io/PGEcore/reference/create_FEM_input.md).

- COI:

  Average complexity of infection.

- threads:

  Number of threads.

- seed:

  Random seed.

- n_chains:

  Number of MCMC chains to run. At least two are needed to compute the
  Gelman-Rubin R-hat convergence diagnostic.

- no_traces_preburnin:

  Number of MCMC traces retained per chain before burn-in is discarded.
  This is the dominant driver of memory use: the sampler pre-allocates a
  `no_traces_preburnin x n_samples x n_haplotypes x n_chains` array of
  doubles, so a 6-locus group (64 haplotypes) at the default reaches
  several GB.

- thinning_interval:

  Number of Metropolis-Hastings updates performed per retained trace.
  Total sampler updates are `no_traces_preburnin * thinning_interval`,
  so halving the traces and doubling this runs the chain exactly as far
  while storing fewer draws. Note that reported ESS is bounded by the
  number of *retained* draws.

## Value

A list with the population frequency table, MCMC runtime, marker names,
alternate alleles, group size, mono-allelic loci, and a
`convergence_diag` data frame.
