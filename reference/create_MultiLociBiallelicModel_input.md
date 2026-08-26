# Create Multi-Loci Biallelic Model input

Create Multi-Loci Biallelic Model input

## Usage

``` r
create_MultiLociBiallelicModel_input(
  input_path,
  loci_group,
  aa_sample_occurence_cut_off = 0
)
```

## Arguments

- input_path:

  Path to amino-acid calls TSV.

- loci_group:

  Path to loci-group TSV (`group_id`, `gene_id`, `aa_position`).

- aa_sample_occurence_cut_off:

  Amino-acid calls must occur in more than this number of samples to be
  included.

## Value

A list (`MLBM_object`) with `MLBM_data`, `staves_data`,
`by_group_table`, and `groups`.
