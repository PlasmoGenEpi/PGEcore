# Translate amino-acid loci from unique haplotypes via overlap alignment

Translate amino-acid loci from unique haplotypes via overlap alignment

## Usage

``` r
translate_microhap_seqs(
  allele_table_unique_haps_tab,
  microhaps_intersected_with_loci_of_interest,
  ref_bed_by_loci_lookup,
  loci_of_interest_tab
)
```

## Arguments

- allele_table_unique_haps_tab:

  Unique `target_name`/`seq` rows.

- microhaps_intersected_with_loci_of_interest:

  Targets covering loci.

- ref_bed_by_loci_lookup:

  Named list of one-row ref_bed tibbles.

- loci_of_interest_tab:

  Loci-of-interest table.

## Value

Tibble of codon/AA calls per haplotype.
