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

```
# Basic usage
Rscript scripts/moire_wrapper/moire_wrapper.R \
    --allele_table data/example2_allele_table.tsv

# Store MCMC results to be able to check convergence
Rscript scripts/moire_wrapper/moire_wrapper.R \
    --allele_table data/example2_allele_table.tsv \
    --mcmc_results_output mcmc_results.rds

# Run with parallel tempering and write per-rung swap acceptance rates
Rscript scripts/moire_wrapper/moire_wrapper.R \
    --allele_table data/example2_allele_table.tsv \
    --num_chains 2 \
    --pt_chains 4 \
    --acceptance_rates_output acceptance_rates.tsv
```

## Outputs

In addition to the five posterior summaries (`--coi_summary`, `--he_summary`,
`--allele_freq_summary`, `--relatedness_summary`, `--effective_coi_summary`), 
the wrapper also writes a convergence diagnostics table (`--convergence_output`) 
and a parallel tempering swap acceptance rates table 
(`--acceptance_rates_output`), if parallel tempering is used.
