# Create a MOIRE input object from an allele table

Reads a TSV of allele presence, validates columns, and packages MCMC
parameters for
[`run_moire()`](https://plasmogenepi.github.io/PGEcore/reference/run_moire.md).

## Usage

``` r
create_moire_input(
  input_path,
  allow_relatedness,
  burnin,
  samples_per_chain,
  thin,
  verbose,
  eps_pos_alpha,
  eps_pos_beta,
  eps_neg_alpha,
  eps_neg_beta,
  r_alpha,
  r_beta,
  mean_coi_shape,
  mean_coi_scale,
  max_eps_pos,
  max_eps_neg,
  record_latent_genotypes,
  pt_chains,
  pt_grad_lower,
  pt_num_threads,
  adapt_temp,
  max_runtime,
  n_chains,
  threads
)
```
