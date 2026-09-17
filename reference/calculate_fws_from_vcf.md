# Calculate within-host Fws from a VCF via moimix

Converts the VCF to GDS with **SeqArray** when the GDS is missing, older
than the VCF, or `overwrite` is `TRUE`, then runs
[`moimix::getFws()`](https://rdrr.io/pkg/moimix/man/getFws.html). The
VCF must carry per-sample allelic depths (`FORMAT/AD`).

## Usage

``` r
calculate_fws_from_vcf(
  vcf,
  output = "fws_result.tsv",
  gds = NULL,
  population_name = NULL,
  overwrite = FALSE,
  verbose = FALSE
)
```

## Arguments

- vcf:

  Input VCF path. See *Inputs*.

- output:

  Output TSV path. Defaults to `"fws_result.tsv"`.

- gds:

  Optional GDS path. If `NULL`, derived from the VCF path.

- population_name:

  Optional population label added as a column.

- overwrite:

  If `TRUE`, rebuild the GDS even when it is up to date.

- verbose:

  If `TRUE`, print progress messages.

## Value

A tibble with `specimen_name`, `fws`, and optionally `population_name`,
sorted by `fws`.

## Details

### Inputs

- **`vcf`**: Input VCF path (`.vcf` or `.vcf.gz`) with `FORMAT/AD`.

### Multi-allelic sites

[`moimix::getFws()`](https://rdrr.io/pkg/moimix/man/getFws.html) is
biallelic-only and does not check: on a record with more than one ALT it
reads only the first two `AD` columns while dividing by the depth over
all of them, and its MAF becomes the rarest allele's frequency. Both
errors bias Fws downward, so multi-allelic input reads as spuriously
polyclonal. Sites with more than two alleles are therefore dropped
before `getFws()` runs, with a warning naming how many. The ALT count is
the one declared in the VCF, so a site whose extra ALT no sample here
has reads for is still dropped; split and left-normalise upstream
(`bcftools norm -m-`) if that loses too many sites. For an Fws that uses
every allele instead of discarding these sites, see
`plasgenomicsutils calculate_fws --multiallelic collapse`.

### Outputs

- **`output`**: Fws TSV with columns `specimen_name`, `fws`, and
  optionally `population_name`. Defaults to `"fws_result.tsv"`.

- **`gds`** (optional): GDS path used/created beside the VCF. If `NULL`,
  derived by replacing `.vcf` / `.vcf.gz` with `.gds`.

### Running

    calculate_fws_from_vcf(
      vcf = "calls.vcf.gz",
      output = "fws_result.tsv"
    )

    Rscript exec/calculate_fws_from_vcf \
      --vcf calls.vcf.gz \
      --output fws_result.tsv

Requires **moimix** and **SeqArray** (Suggests). `moimix` is installed
from GitHub (`bahlolab/moimix`), not CRAN.

## See also

[`vignette("input-formats", package = "PGEcore")`](https://plasmogenepi.github.io/PGEcore/articles/input-formats.md)
