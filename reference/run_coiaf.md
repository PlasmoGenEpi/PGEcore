# Estimate complexity of infection (COI) using coiaf

Processes SNP read-count data and optionally population-level minor
allele frequencies (PLMAF), then estimates COI with both the frequency
and variant methods from the **coiaf** package (Suggests; not installed
with PGEcore).

## Usage

``` r
run_coiaf(snp_calls, plmaf = NULL, seq_error = 0.01, max_coi = 25)
```

## Arguments

- snp_calls:

  SNP-calls data frame. See *Inputs*.

- plmaf:

  Optional PLMAF data frame. See *Inputs*.

- seq_error:

  Sequencing error rate (default: `0.01`).

- max_coi:

  Maximum COI to consider (default: `25`).

## Value

A data frame with columns `specimen_name`, `coi_freq`, and
`coi_variant`.

## Details

### Inputs

- **`snp_calls`**: SNP-calls data frame (`specimen_name`, `snp_name`,
  `reads`, `seq_base`). See
  [`vignette("input-formats", package = "PGEcore")`](https://plasmogenepi.github.io/PGEcore/articles/input-formats.md).

- **`plmaf`**: Optional PLMAF data frame (`snp_name`, `seq_base`,
  `plmaf`). If `NULL`, PLMAF is calculated from `snp_calls`.

### Outputs

- Returns a data frame with `specimen_name`, `coi_freq`, and
  `coi_variant` (not written to disk). For file I/O, use
  [`coiaf_wrapper()`](https://plasmogenepi.github.io/PGEcore/reference/coiaf_wrapper.md).

### Running

    run_coiaf(snp_calls = snp_df, plmaf = plmaf_df)

File and CLI users should call
[`coiaf_wrapper()`](https://plasmogenepi.github.io/PGEcore/reference/coiaf_wrapper.md)
/ `Rscript exec/coiaf_wrapper ...`.

Requires **coiaf** (Suggests).

## See also

[`coiaf_wrapper()`](https://plasmogenepi.github.io/PGEcore/reference/coiaf_wrapper.md),
[`vignette("input-formats", package = "PGEcore")`](https://plasmogenepi.github.io/PGEcore/articles/input-formats.md)
