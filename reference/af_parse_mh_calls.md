# Load microhaplotype calls for naive allele frequency estimation

Load microhaplotype calls for naive allele frequency estimation

## Usage

``` r
af_parse_mh_calls(path)
```

## Arguments

- path:

  Path to a TSV of microhaplotype genotypes.

## Value

A tibble with columns including `specimen_name`, `target_name`, `reads`,
and `variant`.
