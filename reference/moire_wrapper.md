# Run MOIRE from allele-table and output paths

Reads an allele table, runs MOIRE MCMC, and writes COI, He,
allele-frequency, relatedness, effective-COI, and convergence summaries.
Requires **moire**, **checkmate**, and **posterior** (Suggests).

## Usage

``` r
moire_wrapper(
  allele_table,
  allow_relatedness = TRUE,
  burnin = 10000L,
  samples_per_chain = 1000L,
  n_chains = 3L,
  threads = 1L,
  thin = 1L,
  verbose = FALSE,
  eps_pos_alpha = 1,
  eps_pos_beta = 1,
  eps_neg_alpha = 1,
  eps_neg_beta = 1,
  r_alpha = 1,
  r_beta = 1,
  mean_coi_shape = 0.1,
  mean_coi_scale = 10,
  max_eps_pos = 2,
  max_eps_neg = 2,
  record_latent_genotypes = FALSE,
  pt_chains = 1L,
  pt_grad_lower = 0,
  pt_num_threads = 1L,
  adapt_temp = TRUE,
  max_runtime = Inf,
  coi_output = "coi_output.tsv",
  he_output = "he_output.tsv",
  allele_freq_output = "allele_freq_output.tsv",
  relatedness_output = "relatedness_output.tsv",
  effective_coi_output = "effective_coi_output.tsv",
  mcmc_results_output = NULL,
  convergence_output = "convergence_diag.tsv",
  acceptance_rates_output = NULL,
  seed = NULL
)
```

## Arguments

- allele_table:

  Path to allele table TSV. See *Inputs*.

- allow_relatedness:

  Logical; allow relatedness within samples.

- burnin:

  MCMC burn-in iterations.

- samples_per_chain:

  Samples per MCMC chain.

- n_chains:

  Number of independent MCMC chains. More than one is required to
  compute the Gelman-Rubin R-hat convergence diagnostic. This is
  distinct from `pt_chains` (parallel-tempering rungs within a chain).

- threads:

  Threads used to run the independent chains in parallel (MOIRE's
  `num_cores`).

- thin:

  Thinning interval for the MCMC sampler; only every `thin`-th sample is
  retained.

- verbose:

  Logical; verbose MOIRE output.

- eps_pos_alpha, eps_pos_beta:

  Positive error-rate prior.

- eps_neg_alpha, eps_neg_beta:

  Negative error-rate prior.

- r_alpha, r_beta:

  Relatedness prior.

- mean_coi_shape, mean_coi_scale:

  Mean COI prior.

- max_eps_pos, max_eps_neg:

  Maximum error rates.

- record_latent_genotypes:

  Logical; record latent genotypes.

- pt_chains:

  Number of parallel-tempering chains.

- pt_grad_lower:

  Lower bound for PT temperature gradient.

- pt_num_threads:

  Threads for parallel tempering.

- adapt_temp:

  Logical; adaptive temperature.

- max_runtime:

  Maximum MCMC runtime.

- coi_output:

  Output path for COI summary TSV. See *Outputs*.

- he_output:

  Output path for He summary TSV. See *Outputs*.

- allele_freq_output:

  Output path for allele-frequency summary TSV. See *Outputs*.

- relatedness_output:

  Output path for relatedness summary TSV. See *Outputs*.

- effective_coi_output:

  Output path for effective COI summary TSV. See *Outputs*.

- mcmc_results_output:

  Optional RDS path for full MCMC results.

- convergence_output:

  Output path for MCMC convergence diagnostics. See *Outputs*.

- acceptance_rates_output:

  Optional output path for parallel-tempering swap acceptance rates.
  Only meaningful when `pt_chains > 1`.

- seed:

  Integer seed for reproducible sampling, or `NULL` to leave MOIRE
  non-deterministic. Support is detected from
  [`moire::run_mcmc()`](https://EPPIcenter.github.io/moire/reference/run_mcmc.html)'s
  formals; supplying a seed to a MOIRE build that cannot honour it is an
  error rather than a silent fall-through to unseeded sampling.

## Value

Invisibly, the MOIRE MCMC result object.

## Details

### Inputs

- **`allele_table`**: Allele table TSV (`specimen_name`, `target_name`,
  `seq`). See
  [`vignette("input-formats", package = "PGEcore")`](https://plasmogenepi.github.io/PGEcore/articles/input-formats.md).

### Outputs

- **`coi_output`**: COI summary (`specimen_name`, `coi`, …).

- **`he_output`**: Heterozygosity summary (`target_name`, `he`, …).

- **`allele_freq_output`**: Allele frequencies (`target_name`, `seq`,
  `freq`, …).

- **`relatedness_output`**: Within-host relatedness (`specimen_name`,
  `within_host_rel`, …).

- **`effective_coi_output`**: Effective COI (`specimen_name`, `ecoi`,
  …).

- **`convergence_output`**: MCMC diagnostics (`variable`, `mean`,
  `median`, `sd`, `q5`, `q95`, `rhat`, `ess_bulk`, `ess_tail`).

- **`mcmc_results_output`**: Optional RDS of the full MCMC object.

- **`acceptance_rates_output`**: Optional PT swap rates when
  `pt_chains > 1`.

### Running

    moire_wrapper(
      allele_table = "allele_table.tsv",
      coi_output = "coi_output.tsv",
      he_output = "he_output.tsv",
      allele_freq_output = "allele_freq_output.tsv",
      relatedness_output = "relatedness_output.tsv",
      effective_coi_output = "effective_coi_output.tsv"
    )

    Rscript exec/moire_wrapper \
      --allele_table allele_table.tsv \
      --coi_output coi_output.tsv \
      --he_output he_output.tsv \
      --allele_freq_output allele_freq_output.tsv \
      --relatedness_output relatedness_output.tsv \
      --effective_coi_output effective_coi_output.tsv

Requires **moire**, **checkmate**, and **posterior** (Suggests).

## See also

[`run_moire()`](https://plasmogenepi.github.io/PGEcore/reference/run_moire.md),
[`vignette("input-formats", package = "PGEcore")`](https://plasmogenepi.github.io/PGEcore/articles/input-formats.md)
