# Expand STAVE MLAF rows into single-locus allele frequencies

Requires the optional **variantstring** package.

## Usage

``` r
convert_mlaf_to_slaf(dat)
```

## Arguments

- dat:

  MLAF tibble with `group_id`, `variant`, and `freq`.

## Value

Tibble with `group_id`, `gene_id`, `aa_position`, `aa`, and `freq`.
