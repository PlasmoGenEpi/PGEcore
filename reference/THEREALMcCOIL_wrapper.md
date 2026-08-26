# Estimate COI and allele frequencies with THEREALMcCOIL

Runs THEREALMcCOIL MCMC on SNP calls to estimate per-specimen COI and
single-locus allele frequencies. Compiled C routines
(`McCOIL_categorical`, `McCOIL_prop`) are linked at package install
time. Requires **posterior** (Suggests) for convergence diagnostics.

## Usage

``` r
THEREALMcCOIL_wrapper(
  snp_calls,
  slaf_output,
  coi_output,
  model = "categorical",
  maxCOI = 25L,
  threshold_ind = 20L,
  threshold_site = 20L,
  totalrun = 10000L,
  burnin = 1000L,
  M0 = 15L,
  e1 = 0.05,
  e2 = 0.05,
  epsilon = 0.02,
  err_method = 1L,
  seed = 321L,
  n_chains = 3L,
  convergence_output = "convergence_diag.tsv"
)
```

## Arguments

- snp_calls:

  Path to SNP-calls TSV. See *Inputs*.

- slaf_output:

  Output TSV of allele frequencies. See *Outputs*.

- coi_output:

  Output TSV of COI estimates. See *Outputs*.

- model:

  `"categorical"` (heterozygous/homozygous calls) or `"proportional"`
  (allele frequency / read-count data).

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

- convergence_output:

  Output TSV of MCMC convergence diagnostics. See *Outputs*.

## Value

A list with `slaf`, `coi` and `convergence` (invisibly after writing
outputs).

## Details

### Inputs

- **`snp_calls`**: SNP-calls TSV (at least `specimen_name`, `snp_name`,
  `reads`, `seq_base`; typically also `target_name`, `pos`, `he`). See
  [`vignette("input-formats", package = "PGEcore")`](https://plasmogenepi.github.io/PGEcore/articles/input-formats.md).

### Outputs

- **`slaf_output`**: Single-locus allele frequencies (`variant`,
  `freq`).

- **`coi_output`**: COI estimates (`specimen_name`, `coi`).

- **`convergence_output`**: MCMC diagnostics (`variable`, `mean`,
  `median`, `sd`, `q5`, `q95`, `rhat`, `ess_bulk`, `ess_tail`).

### Running

    THEREALMcCOIL_wrapper(
      snp_calls = "snp_calls.tsv",
      slaf_output = "slaf.tsv",
      coi_output = "coi.tsv"
    )

    Rscript exec/THEREALMcCOIL_wrapper \
      --snp_calls snp_calls.tsv \
      --slaf_output slaf.tsv \
      --coi_output coi.tsv

## See also

[`vignette("input-formats", package = "PGEcore")`](https://plasmogenepi.github.io/PGEcore/articles/input-formats.md)
