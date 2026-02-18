# Dcifer Single-Locus Allele Frequency (SLAF) Wrapper

Contents: 
* [Tool Information](#tool-information)
* [Script Usage](#script-usage)

## Tool Information

### Purpose

Dcifer is primarily a method for estimating genetic relatedness between malaria 
samples, but it does contain functionality to estimate the allele frequencies of 
genetic markers with a maximum likelihood method. This is a wrapper for that 
allele frequency functionality. For a wrapper of the genetic relatedness 
estimation, see `scripts/dcifer_ibd_wrapper`.

### Existing resources

Dcifer is described in: Gerlovina, I., Gerlovin, B., Rodríguez-Barraquer, I., & 
Greenhouse, B. (2022). Dcifer: An IBD-based method to calculate genetic distance 
between polyclonal infections. Genetics, 222(2). 
https://doi.org/10.1093/genetics/iyac126

## Input and Output Formats

The script requires an allele table (`--allele_table`) and can either take a 
table with complexity of infection values (`--coi_table`) or can calculate COI 
with `dcifer::getCOI()`. For specifics on input formats, refer to the argument 
documentation in the script.

The output will be a TSV with columns for target\_id, seq, freq, and 
sample\_total.

## Script Usage

```{r}
scripts/dcifer_slaf_wrapper/dcifer_slaf_wrapper.R --allele_table \
    data/example_allele_table.tsv --coi_table data/example_coi_table.tsv \
    --slaf_output slaf.tsv
```
