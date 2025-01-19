# Dcifer

Contents: 
* [Tool Information](#tool-information)
* [Script Usage](#script-usage)

## Tool Information

### Purpose

Dcifer is a method for estimating genetic relatedness between malaria samples 
based on the concept of identity-by-descent (IBD). It is designed for use with 
unlinked loci and can handle polyclonal and multiallelic data. It evaluates the 
significance of each estimate relative to a user-specified null value using a 
likelihood ratio method.

This script wraps the Dcifer tool and produces relatedness estimates for each 
pair of samples in the provided allele table. COI and allele frequency may be 
estimated naively with Dcifer's built-in functions, or provided as input tables.

### Existing resources

Dcifer is described in: Gerlovina, I., Gerlovin, B., Rodríguez-Barraquer, I., & 
Greenhouse, B. (2022). Dcifer: An IBD-based method to calculate genetic distance 
between polyclonal infections. Genetics, 222(2). 
https://doi.org/10.1093/genetics/iyac126

A tutorial and information about this tool can be found 
[here](https://mrc-ide.github.io/PGEforge/tutorials/dcifer/dcifer_background.html).

## Script Usage

To run the wrapper with a user-specified COI table and using two threads:

```{r}
Rscript scripts/dcifer_wrapper/dcifer_wrapper.R --allele_table \
    data/example_allele_table.tsv --coi_table data/example_coi_table.tsv \
    --threads 2 --btwn_host_rel_output btwn_host_rel.tsv
```

To run the wrapper with a user-specified allele frequency table and built-in COI 
estimation:

```{r}
Rscript scripts/dcifer_wrapper/dcifer_wrapper.R --allele_table \
    data/example2_allele_table.tsv --allele_freq_table \
    data/example_slaf_mhap.tsv --btwn_host_rel_output btwn_host_rel.tsv
```
