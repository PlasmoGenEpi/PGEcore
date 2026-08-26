# THEREALMcCOIL proportional MCMC

Ported from `McCOIL_proportional.R`. Reads the fitted beta grid from
package `extdata` and calls compiled `McCOIL_prop`. Named
`run_mccoil_proportional` so it does not collide with registered C
symbols.

## Usage

``` r
run_mccoil_proportional(
  dataA1,
  dataA2,
  maxCOI = 25,
  totalrun = 10000,
  burnin = 1000,
  M0 = 15,
  epsilon = 0.02,
  err_method = 1,
  path = getwd(),
  output = "output.txt"
)
```

## Arguments

- dataA1:

  Allele-1 read counts (samples x sites).

- dataA2:

  Allele-2 read counts (samples x sites).

- maxCOI:

  Upper bound for COI.

- totalrun:

  Total MCMC iterations.

- burnin:

  Burn-in iterations discarded when summarising.

- M0:

  Initial COI.

- epsilon:

  Sequencing error parameter for the proportional model.

- err_method:

  `1`/`2` treat epsilon as constant; `3` estimates it.

- path:

  Directory for MCMC output files.

- output:

  Base filename for MCMC traces.

## Value

`NULL`. Called for the files it writes.
