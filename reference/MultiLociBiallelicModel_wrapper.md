# Estimate multilocus haplotype frequencies with MultiLociBiallelicModel

Estimates multilocus haplotype frequencies from amino-acid calls and
loci groups. Samples with any missing genotype are dropped (legacy
behaviour; the model cannot handle missing data). Requires
**variantstring** 1.x (Suggests).

## Usage

``` r
MultiLociBiallelicModel_wrapper(
  aa_calls,
  loci_groups,
  mlaf_output,
  aa_sample_occurence_cut_off = 0
)
```

## Arguments

- aa_calls:

  Path to amino-acid calls TSV. See *Inputs*.

- loci_groups:

  Path to loci-groups TSV. See *Inputs*.

- mlaf_output:

  Output TSV path. See *Outputs*.

- aa_sample_occurence_cut_off:

  Amino-acid calls must occur in more than this number of samples to be
  included (legacy default `0`).

## Value

The bound MLAF tibble (invisibly after writing `mlaf_output`).

## Details

### Inputs

- **`aa_calls`**: Amino-acid calls TSV (`specimen_name`, `gene_id`,
  `aa_position`, `ref_aa`, `aa`). See
  [`vignette("input-formats", package = "PGEcore")`](https://plasmogenepi.github.io/PGEcore/articles/input-formats.md).

- **`loci_groups`**: Loci-groups TSV (`group_id`, `gene_id`,
  `aa_position`).

### Outputs

- **`mlaf_output`**: Multilocus allele frequencies (`group_id`,
  `variant`, `freq`).

### Running

    MultiLociBiallelicModel_wrapper(
      aa_calls = "aa_calls.tsv",
      loci_groups = "loci_groups.tsv",
      mlaf_output = "mlaf.tsv"
    )

    Rscript exec/MultiLociBiallelicModel_wrapper \
      --aa_calls aa_calls.tsv \
      --loci_groups loci_groups.tsv \
      --mlaf_output mlaf.tsv

Requires **variantstring** (Suggests).

## See also

[`vignette("input-formats", package = "PGEcore")`](https://plasmogenepi.github.io/PGEcore/articles/input-formats.md)
