# Naive Allele Frequency Estimation

Contents: 
* [Tool Information](#tool-information)
* [Script Usage](#script-usage)

## Tool Information

This tool provides a naive implementation of estimating allele frequency from amino acid calls. 

Allele frequency is estimated in two ways:
- Allele frequency by read count: within sample allele proportions are calculated using read counts and then the allele frequency is calculated as the average within sample allele frequency
- Allele frequency by presence/absence: allele frequency is calculated as the count of a particular allele divided by the total number of alleles observed at a given position

## Script Usage 
```
Rscript naive_allele_frequency_wrapper.R \
  --aa_calls data/example_amino_acid_calls.tsv \
  --method read_count_prop \ # or presence_absence
  --output allele_freqs.tsv
```