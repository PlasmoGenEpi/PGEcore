# Naive COI Distribution Estimation

Contents: 
* [Tool Information](#tool-information)
* [Script Usage](#script-usage)

## Tool Information

This tool provides a naive summary of the distribution of COI calls across samples. It takes as input a table of COI calls, with the columns: specimen_id, coi. The output is a table with the number of samples with each COI call and the proportion of samples with that COI call.

## Script Usage 
```
Rscript count_samples_by_coi.R \
  --coi_calls data/example_coi_table.tsv \
  --output coi_distribution.tsv
```