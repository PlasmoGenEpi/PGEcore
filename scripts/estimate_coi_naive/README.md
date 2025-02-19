# Estimate COI via naieve methods

Contents: 
* [Tool Information](#tool-information)
* [Script Usage](#script-usage)

## Tool Information

There are several simple ways to estimate COI from raw data without relying on
models. For example, if a sample contains any loci with 4 alleles then you could
argue that COI is 4 (at least). However, sequencing errors can also create this
false signal. Slightly more nuanced methods that make use of thresholds can be
used to provide results that are robust to errors.

These module implements several basic COI estimation methods. It uses the
following steps:

1. Reads in allele call data, for example produced from a PMO object.
2. For each sample, counts the number of alleles at each locus and sorts these
into decreasing order. Either takes the $n^{th}$ value as the COI estimate (the
"integer method"), or takes the value at a stated quantile (the "quantile
method").
3. Returns a data.frame of COI estimate for each sample.
4. Writes this output data.frame to file.


## Script Usage

The `estimate_coi_naive.R` script contains all the requisite functions to read
in the data, estimate COI, and write results to file.

For example:

```
Rscript estimate_coi_naive.R --input_path data/example2_allele_table.tsv \
    --output_path coi_table.tsv --method integer_method --integer_threshold 5
Rscript estimate_coi_naive.R --input_path data/example2_allele_table.tsv \
    --output_path coi_table.tsv --method quantile_method --quantile_threshold \
    0.05
```
