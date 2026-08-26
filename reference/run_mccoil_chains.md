# Run several independent McCOIL chains into `work_dir`

Prepares the model input once, then runs `n_chains` chains with seeds
`seed, seed + 1, ...`. Chains are independent and write to distinct
trace files in `work_dir`, so they can run in parallel forks; on
Windows, where forking is unavailable, they run sequentially.

## Usage

``` r
run_mccoil_chains(
  df,
  model = "categorical",
  maxCOI = 25,
  threshold_ind = 20,
  threshold_site = 20,
  totalrun = 10000,
  burnin = 1000,
  M0 = 15,
  e1 = 0.05,
  e2 = 0.05,
  epsilon = 0.02,
  err_method = 1,
  seed = 321,
  n_chains = 3,
  work_dir,
  output = "McCOIL_out.txt"
)
```

## Arguments

- df:

  Preprocessed SNP calls.

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

  Random seed for the first chain; chain `i` uses `seed + i - 1`.

- n_chains:

  Number of independent MCMC chains. More than one chain is required for
  the Gelman-Rubin R-hat diagnostic.

- work_dir:

  Directory for McCOIL temp traces (typically under
  [`tempdir()`](https://rdrr.io/r/base/tempfile.html)).

- output:

  Base filename for chain 1, written under `work_dir`; later chains
  append a `_chain<i>` suffix.

## Value

A list with `traces` (per-chain trace data frames, chain 1 first) and
`summary_path`, the summary file written by chain 1.
