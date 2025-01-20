# SNP-Slice README
This package is a wrapper for the main SNP-Slice package, which is available
here:

https://github.com/nianqiaoju/snp-slice/tree/main

Briefly, SNP-Slice takes as inputs a table of reference allele counts and a table
of alternate allele counts (both tables organized as genomic positions in
columns, samples in rows, and read counts at the intersections) and outputs a
table of genotypes (genomic positions in columns, with one row per genotype, and
'1' values for the genomic positions of the input tables that are part of the
genotype of the current row and '0' values for positions that are not part of
the current genotype row). The program also outputs a second table, with one
column for each genotype from the first output table, and one row for each
sample from the input tables, and '1' values for genotypes that are present in a
given sample and '0' values for genotypes that are absent from a sample.

This implementation of SNP-Slice is designed to call SNP-Slice as part of a
nextflow pipeline that analyzes standardized input datasets with a wide variety
of tools. The program first converts the standardized input tables into SNP-Slice
format input tables and
next runs SNP-Slice to generate output tables. Finally, output tables are
reformatted into STAVE format.

## Dependencies
This package currently depends on R and the following R packages:
 - dplyr
 - tidyr
 - readr
 - optparse
 - stringr
 - tibble

## Functions
 - create_SNP_slice_input: converts amino acid tables into SNP-Slice input
 - create_ref_alt: converts standardized input tables into one table per
 amino acid position, with a single reference allele and a single alternate
 allele associated with each position.
 - run_SNPslice_neg_binomial: A function created by Kathryn Murie, borrowed here
 to run the negative binomial model of SNP-Slice.
 - infer_SNPslice_freqs: A second function created by Kathryn Murie to estimate
 the frequencies of genotypes by counting the number of samples that contain
 each genotype divided by the total number of genotype observations in the
 population. Modified slightly to output the number of genotype observations in
 the population.
 - get_alternates: Reads the table created by create_ref_alt into an R data
 frame
 - genotype_to_stave: takes genotypes formatted as 1's and 0's, reference and
 alternate alleles, and original locus amino acid names, and reformats genotypes
 into STAVE format
 - lookup_genotype - looks up whether each locus within the genotype should be
 listed as reference or alternate, as well as the amino acid associated with
 each
 - genotypes_from_freqs - a wrapper function that creates final STAVE formatted
 output datasets