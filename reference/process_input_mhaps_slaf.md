# Read and validate microhaplotype allele frequency tables

Read and validate microhaplotype allele frequency tables

## Usage

``` r
process_input_mhaps_slaf(mhaps_slaf_fnp)
```

## Arguments

- mhaps_slaf_fnp:

  Path to a TSV with columns `target_name`, `seq`, `freq`, and
  `sample_total`.

## Value

A tibble of microhaplotype allele frequencies.
