# Convert a VCF with FORMAT/AD into pileup-style SNP calls

Reverse of
[`snp_calls_to_vcf()`](https://plasmogenepi.github.io/PGEcore/reference/snp_calls_to_vcf.md):
emits one row per `(specimen, SNP, observed allele)` where allelic depth
is at least `min_reads`. Only simple SNP alleles (single-base REF and
ALT in `{A,C,G,T}`) are kept. `pos` is 0-based. `strand` is always `+`.
`target_name` is included only when INFO carries `TARGET=`.

## Usage

``` r
vcf_to_snp_calls(
  vcf,
  snp_calls_output = NULL,
  biallelic = FALSE,
  min_reads = 1L,
  overwrite = FALSE,
  verbose = FALSE
)
```

## Arguments

- vcf:

  Path to a VCF with `FORMAT/AD`. See *Inputs*.

- snp_calls_output:

  Optional output SNP calls TSV path.

- biallelic:

  If `TRUE`, keep only sites with exactly one simple ALT.

- min_reads:

  Minimum AD reads for an allele to be emitted.

- overwrite:

  If `FALSE` (default), refuse to overwrite `snp_calls_output`.

- verbose:

  If `TRUE`, print a summary message when finished.

## Value

A tibble of SNP calls. When `snp_calls_output` is set, the table is also
written to that path.

## Details

### Inputs

- **`vcf`**: Path to a VCF (`.vcf` or `.vcf.gz`) with `FORMAT/AD`.

### Outputs

- **`snp_calls_output`** (optional): SNP calls TSV; gzip-compressed if
  it ends in `.gz`. If `NULL`, results are returned without writing.

### Running

    vcf_to_snp_calls(
      vcf = "calls.vcf.gz",
      snp_calls_output = "snp_calls.tsv"
    )

    Rscript exec/vcf_to_snp_calls \
      --vcf calls.vcf.gz \
      --snp_calls_output snp_calls.tsv

## See also

[`vignette("input-formats", package = "PGEcore")`](https://plasmogenepi.github.io/PGEcore/articles/input-formats.md)
