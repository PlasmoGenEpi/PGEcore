# Convert Multi-Locus Allele Frequency Data to Single-Locus Allele Frequency Data

Contents: 
* [Tool Information](#tool-information)
* [Script Usage](#script-usage)

## Tool Information

Requires: `mrc-ide/variantstring@develop`

This tool provides a converter from multi-locus allele frequency data to single-locus allele frequency data. The input data is a table with columns: group_id, variant, freq. Variant column is in STAVE format. 

The output is a table with columns: group_id, gene, pos, aa, freq.

## Script Usage 
```
Rscript scripts/slaf_from_stave_mlaf/slaf_from_stave_mlaf.R \
  --mlaf_input data/example_mlaf.tsv \
  --output single_locus_allele_frequencies.tsv
```
