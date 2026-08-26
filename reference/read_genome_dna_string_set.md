# Load a FASTA with record names truncated at the first whitespace

Load a FASTA with record names truncated at the first whitespace

## Usage

``` r
read_genome_dna_string_set(path, reason)
```

## Arguments

- path:

  Path to a FASTA file.

- reason:

  Passed to
  [`check_suggested_pkg()`](https://plasmogenepi.github.io/PGEcore/reference/check_suggested_pkg.md)
  for **Biostrings**.

## Value

A `DNAStringSet` with shortened names.
