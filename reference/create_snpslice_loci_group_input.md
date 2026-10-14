# Read loci groups for SNP-Slice multilocus frequencies

Read loci groups for SNP-Slice multilocus frequencies

## Usage

``` r
create_snpslice_loci_group_input(
  loci_groups_path,
  allele_table,
  target_name_col = "target_name",
  allow_multiallelic = FALSE
)
```

## Arguments

- allow_multiallelic:

  Keep targets with more than two alleles. The biallelic observation
  models drop them (SNP-Slice cannot represent a third allele), so by
  default they are excluded from every group with a warning; the
  multinomial model keeps them.
