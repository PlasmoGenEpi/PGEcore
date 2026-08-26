# Extract SNP bases from unique haplotypes via overlap alignment

Extract SNP bases from unique haplotypes via overlap alignment

## Usage

``` r
extract_snps_of_interest(
  allele_table_unique_haps_tab,
  microhaps_intersected_with_snps_of_interest,
  ref_bed_by_loci_lookup,
  snps_of_interest_tab
)
```

## Arguments

- allele_table_unique_haps_tab:

  Unique `target_name`/`seq` rows.

- microhaps_intersected_with_snps_of_interest:

  Targets covering SNPs.

- ref_bed_by_loci_lookup:

  Named list of one-row ref_bed tibbles.

- snps_of_interest_tab:

  SNP-of-interest table.

## Value

Tibble of SNP calls per haplotype.
