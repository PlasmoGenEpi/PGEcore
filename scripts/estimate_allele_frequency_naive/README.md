# Naive Allele Frequency Estimation

Contents: 
* [Tool Information](#tool-information)
* [Script Usage](#script-usage)

## Tool Information

This tool provides a naive implementation of estimating allele frequency from 
either amino acid or microhaplotype calls.

Allele frequency is estimated in two ways:
- Allele frequency by read count: within sample allele proportions are calculated using read counts and then the allele frequency is calculated as the average within sample allele frequency
- Allele frequency by presence/absence: allele frequency is calculated as the count of a particular allele divided by the total number of alleles observed at a given position

## Arguments

`--aa_calls`: A TSV with columns for specimen_id, gene_id, aa_position, aa, and 
read_count.

`--mh_calls`: A TSV with columns for specimen_id, target_id, seq, and 
read_count.

`--method`: A string containing either "read_count_prop" or "presence_absence".

`--output`: Path for output TSV. If `--aa_calls` was provided, it will have 
columns for variant, formatted as a STAVE string, and freq. If `--mh_calls` was 
provided, it will have columns for target_id, seq, and freq.

## Example Usage 

```
# Example with amino acid calls and read count proportion method
Rscript scripts/estimate_allele_frequency_naive/estimate_allele_frequency_naive.R \
  --aa_calls data/example_amino_acid_calls.tsv \
  --method read_count_prop \
  --output allele_freqs.tsv
# Example with microhaplotypes and the presence/absence method
Rscript scripts/estimate_allele_frequency_naive/estimate_allele_frequency_naive.R \
  --aa_calls data/example_allele_table.tsv \
  --method presence_absence \
  --output allele_freqs.tsv
```
