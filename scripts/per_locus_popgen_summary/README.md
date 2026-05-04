# Allele per Locus Summary

Contents:

* [Tool Information](#tool-information)
* [Script Usage](#script-usage)

## Tool Information

This tool provides a simple calculation of popgen stats. By locus, it calculates:

* Number of Segregating Sites.
* Nucleotide Diversity.
* Tajima's D.

## Script Usage 

```
Rscript scripts/per_locus_popgen_summary/per_locus_popgen_summary.R \
    --allele_table <allele_table.tsv>
```
