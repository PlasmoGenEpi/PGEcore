# Population-genetic statistics for one locus

Aligns unique allele sequences with **msa**, restores duplicate
sequences, then computes nucleotide diversity, segregating sites, and
Tajima's D with **pegas**. Identical sequences skip alignment and return
zeros.

## Usage

``` r
calculate_popgen_stats(allele_data, msa_method = "Muscle")
```

## Arguments

- allele_data:

  Character vector of allele sequences.

- msa_method:

  Passed to [`msa::msa()`](https://rdrr.io/pkg/msa/man/msa.html):
  `"Muscle"`, `"ClustalW"`, or `"ClustalOmega"`.

## Value

A named list of statistics.

## Details

Requires the **msa** alignment binary on `PATH` (`muscle`, `clustalw`,
or `clustalo` depending on `msa_method`). The **msa** R package does not
bundle those tools; install them separately (for example via Conda).
