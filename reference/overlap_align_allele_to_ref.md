# Overlap-align an allele DNAString to a reference sequence

Overlap-align an allele DNAString to a reference sequence

## Usage

``` r
overlap_align_allele_to_ref(allele_seq, ref_seq, mat)
```

## Arguments

- allele_seq:

  A `DNAString` (already oriented to the plus strand of the target).

- ref_seq:

  Character reference sequence (already oriented).

- mat:

  Substitution matrix from
  [`nucleotide_overlap_substitution_matrix()`](https://plasmogenepi.github.io/PGEcore/reference/nucleotide_overlap_substitution_matrix.md).

## Value

A pairwise alignment object.
