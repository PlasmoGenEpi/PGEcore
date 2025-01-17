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
  are present. Doing so would require running more complex methods that attempt to phase genotypes.

We can use these rules to obtain all known phased genotypes from the raw data.
For example, imagine we have the following two samples, defined in
[variant string format](https://github.com/mrc-ide/variantstring):

```
crt:1_2:A_A/C
crt:1_2_3:A/C_A_A/C
```

The first sample we can unambiguously phase into the following two genotypes:

```
crt:1_2:A_A
crt:1_2:A_C
```

The second sample cannot be unambiguously phased due to the two heterozygous loci.

Now imagine we are comparing both of these phased genotypes back against the
original two samples. We would obtain the following matches:

```
crt:1_2:A_A --> present unambiguously in sample1 and sample2
crt:1_2:A_C --> present unambiguously in sample1 only
```

Notice that for sample2, even though we were not able to extract its component
genotypes, we were still able to identify unambiguous matches against this
sample. This is because we were not looking over all three positions, we were
only looking at two positions in which there was only one heterozygous locus.

This script performs these operations over an entire dataset. Briefly, it takes
the following steps:

1. Convert input data into variant string format
2. Extract all unambiguous phased genotypes from the data
3. Compare all phased genotypes against the data to estimate prevalence and frequency

Genotype frequency in this case is obtained as the average of the within-sample
genotype frequencies (WSGF) over all samples. This has the advantage of being an
unbiased estimator of the true population-level genotype frequency (PLGF)
without requiring a known complexity of infection for each sample.


## Script Usage

The `multilocus_prevfreq_naive.R` script contains all the requisite functions to
read in the data, calculate prevalence and frequency, and write results to file.

An example of usage, executed from the root of this repo, would be:

``
Rscript scripts/multilocus_prevfreq_naive/multilocus_prevfreq_naive.R \
    --input_path data/example_amino_acid_calls.tsv --output_path mlafp.tsv
``
