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
  M0 = "1,5,15",
  e1 = 0.05,
  e2 = 0.05,
  err_method = 1L,
  seed = 321L,
  n_chains = 3L,
  threads = 1L,
  convergence_output = "convergence_diag.tsv",
  convergence_summary_output = "convergence_summary.tsv"
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
  (allele read counts, with each sample's strain proportions shared
  across loci and read noise estimated).

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

  Starting COI: one value for every chain, or one per chain as a vector
  or comma-separated string (default `"1,5,15"`, for 3 chains). Spread
  starts let R-hat reveal a chain that has not converged.

- e1:

  Probability of calling homozygous loci heterozygous (categorical).

- e2:

  Probability of calling heterozygous loci homozygous (categorical).

- err_method:

  `1`: treat `e1`/`e2` as constants; `3`: estimate them with COI and
  allele frequencies (categorical model).

- seed:

  Random seed for the first chain; chain `i` uses `seed + i - 1`.

- n_chains:

  Number of independent MCMC chains. More than one chain is required for
  the Gelman-Rubin R-hat diagnostic.

- threads:

  Cores used to run chains simultaneously (capped at `n_chains`).

- convergence_output:

  Output TSV of MCMC convergence diagnostics. See *Outputs*.

- convergence_summary_output:

  Output TSV of the run-level convergence summary. See *Outputs*.

## Value

A list with `slaf`, `coi` and `convergence` (invisibly after writing
outputs).

## Details

### Inputs

- **`snp_calls`**: SNP-calls TSV (at least `specimen_name`, `snp_name`,
  `reads`, `ref_base`, `seq_base`; typically also `target_name`, `pos`,
  `he`). See
  [`vignette("input-formats", package = "PGEcore")`](https://plasmogenepi.github.io/PGEcore/articles/input-formats.md).

### Outputs

- **`slaf_output`**: Single-locus allele frequencies (`variant`,
  `snp_name`, `ref_base`, `seq_base`, `freq`). One row per SNP giving
  the frequency of its non-reference allele, with `variant` written as
  `snp_name:seq_base`. The reference allele's frequency is `1 - freq`.
  Both models use only biallelic SNPs whose reference base is observed;
  other SNPs are dropped with a warning.

- **`coi_output`**: COI estimates (`specimen_name`, `coi`).

Both are posterior medians over the post-burn-in draws of all chains.

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
