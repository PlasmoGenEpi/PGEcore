# Estimate multi-locus prevalence and frequency via naive methods

Contents: 

* [Tool Information](#tool-information)
* [Script Usage](#script-usage)

## Tool Information

The prevalence of a multi-locus genotype is defined as the proportion of samples
in which it is detected. In contrast, genotype frequency, while less
straightforward to define, represents the probability that a new malaria
infection carries this specific genotype. Prevalence and frequency are not
equivalent because individuals can harbour multiple strains of malaria
simultaneously. As a result, the sum of prevalence values across all genotypes
may exceed 1, whereas genotype frequencies must always sum to exactly 1.

Calculating genotype prevalence and frequency is more complex than it first
appears, and depends on the number of heterozygous loci:

- If there are **zero heterozygous loci** then we know exactly which genotype is present.
  The genotype is fully phased.
- If there is **a single heterozygous locus** then we know which two phased genotypes must
  be present in the sample. From the relative read counts at this heterozygous locus we can obtain a
  rough estimate of the within-sample proportions of each genotype.
- If there are **two or more heterozygous loci** then we cannot unambiguously state which genotypes
  are present. Doing so would require running more advanced methods that attempt to phase genotypes. However these loci can be filtered with a reasonably high within sample frequency to infer a phased genotype despite the multiple heterogyous sites. E.g. if all loci have an allele at >=70% then it is highly probable that a phased genotype containing all those alleles is present. 

We can use these rules to obtain all known phased genotypes from the raw data.
For example, imagine we have the following two samples, defined in
[variant string format](https://github.com/mrc-ide/variantstring):

This script performs these operations over a set of sample for a given group of loci. Briefly, it takes the following steps:

1. Filter input loci data to the loci of interest in each group
2. Determine unambiguous phased genotypes from the data by: 
	*  Determining samples with only a single heterozygous locus, and combining all variants at that site  
	*  Filtering all loci to a within sample frequency (default is 0.70 and can be changed with \-\-wsaf\_cut\_off) and if all loci have 1 called allele then adding this as a phased multilocus haplotype  


## Script Usage

The `multilocus_prevfreq_naive.R` script contains all the functions needed to
read in the data, calculate prevalence and frequency, and write results to file.

An example of usage, executed from the root of this repo, would be:

```
Rscript ${projectDir}/bin/PGEcore/scripts/multilocus_prevfreq_naive/multilocus_prevfreq_naive.R \
        --aa_table data/example_amino_acid_calls.tsv \
        --loci_groups_input data/example_loci_groups.tsv\
        --output_path mlafp.tsv 

# You can also export single locus allele frequency and prevalences re-calculated from the multilocus calls, this can be useful for comparing to these measures calculated directly off the data for a sanity check of the multilocus calls e.g. if they are extremely different than the calculations made directly from the data than the multilocus calls may be missing important haplotypes 

Rscript ${projectDir}/bin/PGEcore/scripts/multilocus_prevfreq_naive/multilocus_prevfreq_naive.R \
        --aa_table data/example_amino_acid_calls.tsv \
        --loci_groups_input data/example_loci_groups.tsv\
        --output_path mlafp.tsv \
        --recalc_single_locus_output_path aa_sl_from_ml.tsv"
```
