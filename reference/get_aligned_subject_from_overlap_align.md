# Aligned subject sequence including end gaps from an overlap alignment

Aligned subject sequence including end gaps from an overlap alignment

## Usage

``` r
get_aligned_subject_from_overlap_align(pw_overlapAlign)
```

## Arguments

- pw_overlapAlign:

  A
  [`pwalign::pairwiseAlignment()`](https://rdrr.io/pkg/pwalign/man/pairwiseAlignment.html)
  overlap result.

## Value

Character string of the subject with terminal gaps relative to the
pattern.
