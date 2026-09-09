# Lin's concordance correlation coefficient

Computes Lin's (1989) concordance correlation coefficient (CCC) for
agreement between two sets of paired measurements. Unlike Pearson's
correlation, the CCC combines precision (tightness of the points about
their best-fit line) and accuracy (how far that line deviates from the
45-degree line of perfect concordance), so it measures reproducibility
rather than linear association. Values range from -1 to 1, with 1
indicating perfect agreement. Used here to compare per-host COI
estimates between SNP-Slice restarts.

## Usage

``` r
snpslice_ccc(x, y)
```

## Arguments

- x:

  Numeric vector, the first set of measurements.

- y:

  Numeric vector, the second set of measurements, same length as `x`.

## Value

A single numeric value, the concordance correlation coefficient.

## Details

Only pairs where both `x` and `y` are non-missing are used. Returns
`NA_real_` for fewer than three complete pairs and 1 when both vectors
are constant and equal.

## References

Lin L (1989). A concordance correlation coefficient to evaluate
reproducibility. *Biometrics* 45: 255-268.

Lin L (2000). A note on the concordance correlation coefficient.
*Biometrics* 56: 324-325.
