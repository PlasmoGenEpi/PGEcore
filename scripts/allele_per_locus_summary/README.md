# Allele per Locus Summary

Contents: 
* [Tool Information](#tool-information)
* [Script Usage](#script-usage)

## Tool Information

This tool provides a simple count of alleles in the dataset. By locus, it counts:
* Total Allele Count: The total number of alleles (length(allele)).
* Unique Allele Count: The count of unique alleles (length(unique(allele))).
* Allele Singlets: The number of alleles that appear only once sum(table(allele) == 1)  

## Script Usage 
```
Rscript allele_per_locus_summary.R --allele_table <allele_table.tsv>
```
