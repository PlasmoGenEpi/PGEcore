# Naive Allele Frequency and Prevalence Estimation

Contents: 
* [Tool Information](#tool-information)
* [Script Usage](#script-usage)

## Tool Information

This tool provides a naive implementation of estimating allele prevalence from amino acid calls. Prevalence is estimated by calculating the proportion of samples (with at least one call) that have a given variant at a given position. 

## Script Usage 
```
Rscript estimate_allele_prevalence_naive.R \
  --aa_calls data/example_amino_acid_calls.tsv \
  --output prevalence.tsv
```