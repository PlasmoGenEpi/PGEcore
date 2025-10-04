# Incompelete data model (IDM)

Contents:

- [Tool Information](#tool-information)
- [Script Usage](#script-usage)

## Tool Information

| Tool Summary    |                                                                                                                                                                                                                 |
| --------------- | --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| Main use-cases  | Estimate allele frequency and MOI accounting for incomplete data                                                                                                                                                |
| Authors         | Meraj Hashemi, Kristan Schneider                                                                                                                                                                                |
| Latest version  | Unknown                                                                                                                                                                                                         |
| License         | Unknown                                                                                                                                                                                                         |
| Website         | [script](https://doi.org/10.1371/journal.pone.0287161.s002) and [documentation](https://doi.org/10.1371/journal.pone.0287161.s003)                                                                              |
| Code repository | https://github.com/Maths-against-Malaria/MOI---Incomplete-Data-Model.git (not active)                                                                                                                           |
| Publication     | Hashemi M, Schneider KA (2024) Estimating multiplicity of infection, allele frequencies, and prevalences accounting for incomplete data. PLoS ONE 19(3): e0287161. https://doi.org/10.1371/journal.pone.0287161 |

### Purpose

This tool provides a statistical model to estimate MOI and lineage
frequencies/prevalences from single-locus molecular data characterized by
incomplete information.

A note from the [original paper](https://doi.org/10.1371/journal.pone.0287161):

> the new method is recommendable only for data sets in which the molecular
> assays produced poor-quality results. This will be particularly true if the
> model is extended to accommodate information from multiple molecular markers
> at the same time, and incomplete information at one or more markers leads to a
> strong depletion of sample size.

### Existing Resources

- The existing tutorial is in PDF format as mentioned above
  (https://doi.org/10.1371/journal.pone.0287161.s003) and is from in the [original
  paper](https://doi.org/10.1371/journal.pone.0287161). If you find another one,
  please contribute by filing an issue or submitting a pull request.

- Related paper:  
   Hashemi M, Schneider KA. Bias-corrected maximum-likelihood estimation of
  multiplicity of infection and lineage frequencies. PLOS ONE. 2022; 16(12):1–28.
  (This paper describes the original model (OM))

## Script Usage

Assuming you are in the directory containing the wrapper script
(`IDM_wrapper.R`), you can run the script as follows.

```sh
# use allele table as input
 Rscript ./IDM_wrapper.R \
  --allele_table_input ../../data/example_allele_table.tsv \
  --model IDM \
  --slaf_output out_slaf.tsv \
  --eps_initial 0.1 \
  --lambda_initial 0.1
# use amino acide calls as input
 Rscript ./IDM_wrapper.R \
  --aa_calls_input ../../data/example_amino_acid_calls.tsv \
  --model IDM \
  --slaf_output out_slaf.tsv \
  --eps_initial 0.1 \
  --lambda_initial 0.1
```

Note that in this wrapper code,

- Input data is an allele table (`--allele_table_input`, see
  `../../data/example_allele_table.tsv` for example), or
   a table of amino acid calls. See 
  `../../data/example_amino_acid_calls.tsv` for example.
- Output data are presented in a table of allele frequencies. Only lineage (allele) frequencies are reported, while other
  estimates like MOI and prevalence counts are not included. Command line options
  are available to adjust parameters for the maximum-likelihood estimator from the
  tool, including the `--model [OM|IDM]` option, where OM represents the original
  model and IDM represents the incomplete-data model.
