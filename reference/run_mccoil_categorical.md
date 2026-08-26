# THEREALMcCOIL categorical MCMC

Ported from `McCOIL_categorical.R`. Calls compiled `McCOIL_categorical`
in the PGEcore shared library. Writes MCMC traces and a summary TSV
under `path`. Named `run_mccoil_categorical` so it does not collide with
the registered C entry point `McCOIL_categorical`.

## Usage

``` r
run_mccoil_categorical(
  data,
  maxCOI = 25,
  threshold_ind = 20,
  threshold_site = 20,
  totalrun = 10000,
  burnin = 1000,
  M0 = 15,
  e1 = 0.05,
  e2 = 0.05,
  err_method = 1,
  path = getwd(),
  output = "output.txt"
)
```

## Arguments

- data:

  Numeric matrix of heterozygous/homozygous scores (samples x sites);
  missing coded as `-1`.

- maxCOI:

  Upper bound for COI.

- threshold_ind:

  Minimum non-missing sites for a sample.

- threshold_site:

  Minimum non-missing samples for a site.

- totalrun:

  Total MCMC iterations.

- burnin:

  Burn-in iterations discarded when summarising.

- M0:

  Initial COI.

- e1:

  Probability of calling homozygous loci heterozygous.

- e2:

  Probability of calling heterozygous loci homozygous.

- err_method:

  `1`/`2` treat error rates as constants; `3` estimates them.

- path:

  Directory for MCMC output files.

- output:

  Base filename for MCMC traces (`<output>_summary.txt` is also
  written).

## Value

`NULL`. Called for the files it writes.
