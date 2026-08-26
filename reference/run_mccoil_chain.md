# Run a single McCOIL chain and return its per-iteration trace

Runs one MCMC chain of the requested model with a given seed, writing
its trace and summary to `output` under `work_dir`, and reads the trace
back.

## Usage

``` r
run_mccoil_chain(
  mccoil_input,
  model,
  maxCOI,
  threshold_ind,
  threshold_site,
  totalrun,
  burnin,
  M0,
  e1,
  e2,
  epsilon,
  err_method,
  seed,
  work_dir,
  output
)
```

## Arguments

- mccoil_input:

  Prepared model input: the genotype matrix for the categorical model,
  or a list with `a1`/`a2` matrices for the proportional model.

- model:

  `"categorical"` or `"proportional"`.

- maxCOI:

  Upper bound for COI.

- threshold_ind:

  Minimum sites per sample (categorical model).

- threshold_site:

  Minimum samples per locus (categorical model).

- totalrun:

  Total MCMC iterations.

- burnin:

  Burn-in iterations.

- M0:

  Initial COI.

- e1:

  Probability of calling homozygous loci heterozygous (categorical).

- e2:

  Probability of calling heterozygous loci homozygous (categorical).

- epsilon:

  Error parameter for the proportional model.

- err_method:

  `1`: treat error rates as constants; `3`: estimate them with COI and
  allele frequencies.

- seed:

  Random seed for this chain.

- work_dir:

  Directory for McCOIL temp traces.

- output:

  Trace filename (under `work_dir`) for this chain.

## Value

A data frame of the raw per-iteration trace, one row per iteration plus
the trailing acceptance-count row.
