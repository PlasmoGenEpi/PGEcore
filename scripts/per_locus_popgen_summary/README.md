# Allele per Locus Summary

Contents:

* [Tool Information](#tool-information)
* [Script Usage](#script-usage)

## Tool Information

This tool provides a simple calculation of popgen stats. By locus, it calculates:

* Nucleotide diversity (`nucleotide_diversity`)
* Number of segregating sites (`segregating_sites`)
* Tajima's D (`tajima_d`), with *p*-values calculated according to both a normal 
  distribution (`tajima_d_pval_normal`) and a beta distribution 
  (`tajima_d_pval_beta`). See `pegas::tajima.test()` for details.

## Script Usage 

```
Rscript scripts/per_locus_popgen_summary/per_locus_popgen_summary.R \
    --allele_table data/example2_allele_table.tsv \
    --out popgen_summary.tsv
```
