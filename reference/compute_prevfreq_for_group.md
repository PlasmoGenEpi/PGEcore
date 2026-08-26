# Compute prev/freq for one loci group via variantstring

Compute prev/freq for one loci group via variantstring

## Usage

``` r
compute_prevfreq_for_group(group_loci, aa_table)
```

## Arguments

- group_loci:

  Tibble with `gene_id` and `aa_position`.

- aa_table:

  Amino acid calls in long form.

## Value

Tibble with `variant`, `prev`, `freq`, `sample_total`.
