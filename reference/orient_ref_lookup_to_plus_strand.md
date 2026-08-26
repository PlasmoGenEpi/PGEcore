# Reverse-complement `ref_seq` for minus-strand targets in a lookup list

Operates on a copy so repeated calls do not keep reverse-complementing.

## Usage

``` r
orient_ref_lookup_to_plus_strand(ref_bed_by_loci_lookup)
```

## Arguments

- ref_bed_by_loci_lookup:

  Named list of one-row ref_bed tibbles.

## Value

A copy of the lookup with minus-strand `ref_seq` reverse-complemented.
