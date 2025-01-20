# Moire

Contents: 
* [Tool Information](#tool-information)
* [Script Usage](#script-usage)

## Tool Information

A tutorial and information about this tool can be found [here](https://mrc-ide.github.io/PGEforge/tutorials/moire/moire_background.html)

### Purpose

The moire (Multiplicity Of Infection and allele frequency REcovery) tool can be used to estimate allele frequencies, MOI, and within-host relatedness from genetic data subject to experimental error. It utilises a Markov Chain Monte Carlo (MCMC) based approach to Bayesian estimation and can take both polyallelic and SNP data as inputs. This tool also introduces a new metric called effective MOI (eMOI), which combines MOI and within-host relatedness into a unified and comparable measure of genetic diversity.

### Existing resources

The [moire website](https://eppicenter.github.io/moire/index.html) provides basic usage instructions.
Within the [moire website](https://eppicenter.github.io/moire/articles/mcmc_demo.html) there is a more in depth tutorial using simulated genotyping data.

## Script Usage

The script can be run using the following command and default parameters

```
Rscript moire_wrapper.R --allele_table <allele_table.tsv>
```

Running the with all available parameters:

```
Rscript Moire_wrapper.R --allele_table <allele_table.tsv> \
  --allow_relatedness TRUE \
  --burnin 10000 \
  --samples_per_chain 1000 \
  --verbose FALSE \
  --eps_pos_alpha 1 \
  --eps_pos_beta 1 \
  --eps_neg_alpha 1 \
  --eps_neg_beta 1 \
  --r_alpha 1 \
  --r_beta 1 \
  --mean_coi_shape 0.1 \
  --mean_coi_scale 10 \
  --max_eps_pos 2 \
  --max_eps_neg 2 \
  --record_latent_genotypes FALSE \
  --num_chains 1 \
  --num_cores 1 \
  --pt_chains 1 \
  --pt_grad 1 \
  --pt_num_threads 1 \
  --adapt_temp FALSE \
  --max_runtime Inf \
  --seed 1 \
  --coi_summary <coi_summary.tsv> \ 
  --he_summary <he_summary.tsv> \
  --allele_freq_summary <allele_freq_summary.tsv> \
  --relatedness_summary <relatedness_summary.tsv> \
  --effective_coi_summary <effective_coi_summary.tsv> 
```

