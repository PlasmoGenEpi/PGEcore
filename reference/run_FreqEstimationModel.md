# Run FreqEstimationModel MCMC for one group

Run FreqEstimationModel MCMC for one group

## Usage

``` r
run_FreqEstimationModel(sample_matrix_list, COI, threads, seed, n_chains = 3L)
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

## Value

A list with the population frequency table, MCMC runtime, marker names,
alternate alleles, group size, mono-allelic loci, and a
`convergence_diag` data frame.
