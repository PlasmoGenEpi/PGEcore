# Build naive multilocus haplotype calls for one loci group

Build naive multilocus haplotype calls for one loci group

## Usage

``` r
build_naive_multilocus_calls_for_group(aa_calls, group_df, wsaf_cut_off)
```

## Arguments

- aa_calls:

  Amino acid calls with `n_aa`.

- group_df:

  Loci-group rows including `loci_in_group`.

- wsaf_cut_off:

  Dominant-allele WSAF threshold.

## Value

A tibble of specimen haplotypes (`specimen_name`, `variant`, `wsaf`).
