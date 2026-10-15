# THEREALMcCOIL proportional MCMC

Calls compiled `McCOIL_prop` and writes its trace and a
`<output>_summary.txt` in the same layout as
[`run_mccoil_categorical()`](https://plasmogenepi.github.io/PGEcore/reference/run_mccoil_categorical.md),
with `rho` and `pout` as extra parameters.

## Usage

``` r
run_mccoil_proportional(
  dataA1,
  dataA2,
  maxCOI = 25,
  totalrun = 10000,
  burnin = 1000,
  M0 = 1,
  err = 1e-04,
  path = getwd(),
  output = "output.txt"
)
```

## Arguments

- dataA1:

  Allele-1 read counts (samples x sites); negative marks a missing
  locus.

- dataA2:

  Allele-2 read counts (samples x sites).

- maxCOI:

  Upper bound for COI.

- totalrun:

  Total MCMC iterations.

- burnin:

  Burn-in iterations discarded when summarising.

- M0:

  Starting COI for every sample.

- err:

  Per-read error rate.

- path:

  Directory for MCMC output files.

- output:

  Base filename for the MCMC trace.

## Value

`NULL`. Called for the files it writes.

## Details

Each sample has one set of strain proportions shared by all its loci,
and each strain carries one allele per locus. Read noise (`rho`) and the
outlier weight (`pout`) are estimated. A minor allele below the
detection limit, taken as the smallest minor-allele fraction and read
count seen in the data, is treated as censored.
