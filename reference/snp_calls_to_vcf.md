# Build a VCF from pileup SNP calls

Converts pileup SNP calls (raw or collapsed) into a VCF with per-sample
allelic depths (`FORMAT/AD`). Reads are summed per
`(specimen, SNP, allele)` across overlapping targets. `REF`/`ALT` are
written on the forward (genome) strand. Monomorphic sites are skipped.

## Usage

``` r
snp_calls_to_vcf(
  snp_calls,
  genome,
  vcf_output,
  biallelic = FALSE,
  ploidy = 2L,
  gt_min_reads = 1L,
  overwrite = FALSE,
  verbose = FALSE
)
```

## Arguments

- snp_calls:

  Path to a SNP calls TSV, or a data frame with the same columns. See
  *Inputs*.

- genome:

  Path to a reference genome FASTA. See *Inputs*.

- vcf_output:

  Output VCF path; gzip-compressed if it ends in `.gz`.

- biallelic:

  If `TRUE`, keep only sites with a single ALT.

- ploidy:

  Ploidy used to render the GT field.

- gt_min_reads:

  Minimum reads for an allele to count as present in GT.

- overwrite:

  If `FALSE` (default), refuse to overwrite `vcf_output`.

- verbose:

  If `TRUE`, print a summary message when finished.

## Value

Invisibly returns `TRUE`.

## Details

### Inputs

- **`snp_calls`**: SNP calls (`specimen_name`, `chrom`, `pos`,
  `snp_name`, `strand`, `ref_base`, `seq_base`, `reads`), as a file path
  or data frame. `pos` is 0-based (emitted as 1-based VCF `POS`). See
  [`vignette("input-formats", package = "PGEcore")`](https://plasmogenepi.github.io/PGEcore/articles/input-formats.md).

- **`genome`**: Reference genome FASTA used for `##contig` lengths.

### Outputs

- **`vcf_output`**: Output VCF path; gzip-compressed if it ends in
  `.gz`.

### Running

    snp_calls_to_vcf(
      snp_calls = "snp_calls.tsv",
      genome = "genome.fasta",
      vcf_output = "calls.vcf.gz"
    )

    Rscript exec/snp_calls_to_vcf \
      --snp_calls snp_calls.tsv \
      --genome genome.fasta \
      --vcf_output calls.vcf.gz

Requires **Biostrings** (Suggests) to read `--genome` contig lengths.

## See also

[`vignette("input-formats", package = "PGEcore")`](https://plasmogenepi.github.io/PGEcore/articles/input-formats.md)
