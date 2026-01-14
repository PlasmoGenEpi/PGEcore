# Naive Allele Prevalence Estimation

Contents: 
* [Tool Information](#tool-information)
* [Script Usage](#script-usage)

## Tool Information

This tool provides a naive implementation of estimating allele prevalence from 
amino acid calls. Prevalence is estimated by calculating the proportion of 
samples (with at least one call) that have a given variant at a given position. 
Two input formats are accepted: a table of amino acid calls or a table of 
microhaplotype calls.

## Arguments

`--aa_calls`: A TSV with columns for specimen_id, gene_id, aa_position, and aa.

`--mh_calls`: A TSV with columns for specimen_id, target_id, and seq.

`--output`: Path for output TSV. If `--aa_calls` was provided, it will have 
columns for variant, formatted as a STAVE string, prev, and sample_total. If 
`--mh_calls` was provided, it will have columns for target_id, seq, prev, and 
sample_total.

## Script Usage 
```
# Example with amino acid calls
estimate_allele_prevalence_naive.R \
  --aa_calls data/example_amino_acid_calls.tsv \
  --output prevalence.tsv
# Example with microhaplotypes
estimate_allele_prevalence_naive.R \
  --mh_calls data/example_allele_table.tsv \
  --output prevalence.tsv
```
