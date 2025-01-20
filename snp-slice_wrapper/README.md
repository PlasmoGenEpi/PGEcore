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
 - alfy_create_ref_alt: converts standardized input tables into separate
 alternate and reference counts.
 - alfy_create_SNP_slice_input: converts amino acid tables into SNP-Slice input
 - run_SNPslice_neg_binomial
 - infer_SNPslice_freqs
 - get_alternates
 - genotype_to_stave
 - lookup_genotype
 - genotypes_from_freqs