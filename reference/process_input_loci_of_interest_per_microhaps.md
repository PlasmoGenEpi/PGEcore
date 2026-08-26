# Read and validate translated loci-of-interest for microhaplotypes

Read and validate translated loci-of-interest for microhaplotypes

## Usage

``` r
process_input_loci_of_interest_per_microhaps(
  loci_of_interest_per_microhaps_fnp
)
```

## Arguments

- loci_of_interest_per_microhaps_fnp:

  Path to a TSV with columns `target_name`, `gene_id`, `aa_position`,
  `seq`, and `aa`.

## Value

A tibble of translated loci.
