# Derive a fixed-ploidy GT string from per-allele depths

Alleles with depth `>= min_reads` are "present". They fill the `ploidy`
slots ordered by depth: more present than ploidy keeps the top `ploidy`;
fewer pads with the most-supported present allele; none present yields
all-missing (`./.`). Returned indices are 0-based (`0` = REF) and
sorted, joined by `/`.

## Usage

``` r
derive_gt(depths, ploidy, min_reads)
```

## Arguments

- depths:

  Integer vector of depths in REF, ALT order.

- ploidy:

  Ploidy used to render GT.

- min_reads:

  Minimum reads for an allele to count as present.

## Value

A GT string such as `"0/1"` or `"./."`.
