# Collapse identical dictionary rows of a SNP-Slice allocation

In original implementation by Ju et al. (2024), identical SNP haplotypes
are counted once when computing COI, allele frequencies and
heterozygosity: duplicated dictionary rows are removed before any
estimate is formed. SNP-Slice can hold two active strains with the same
haplotype, and a specimen assigned both would then contribute two to its
COI and two strain copies to every frequency. This merges such strains:
a specimen carries the merged strain if it carried any of the
duplicates.

## Usage

``` r
snpslice_dedup_matrices(A, D)
```

## Arguments

- A:

  Allocation matrix, specimens x strains.

- D:

  Dictionary matrix, strains x targets.

## Value

List with the merged `A` and `D`; strains keep first-appearance order.
