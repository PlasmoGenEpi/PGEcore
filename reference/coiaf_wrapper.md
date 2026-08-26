# Estimate COI with coiaf from SNP-call and output paths

Reads SNP calls (and optional PLMAF), runs
[`run_coiaf()`](https://plasmogenepi.github.io/PGEcore/reference/run_coiaf.md),
and writes COI estimates. Requires **coiaf** (Suggests).

## Usage

``` r
coiaf_wrapper(snp_calls, output, plmaf = NULL, seq_error = 0.01, max_coi = 25)
```

## Arguments

- snp_calls:

  Path to SNP-calls TSV. See *Inputs*.

- output:

  Output TSV path. See *Outputs*.

- plmaf:

  Optional path to PLMAF TSV. See *Inputs*.

- seq_error:

  Sequencing error rate (default: `0.01`).

- max_coi:

  Maximum COI to consider (default: `25`).

## Value

The result tibble (also written to `output`).

## Details

### Inputs

- **`snp_calls`**: Path to SNP-calls TSV. See
  [`vignette("input-formats", package = "PGEcore")`](https://plasmogenepi.github.io/PGEcore/articles/input-formats.md).

- **`plmaf`**: Optional path to PLMAF TSV (`snp_name`, `seq_base`,
  `plmaf`).

### Outputs

- **`output`**: TSV with `specimen_name`, `coi_freq`, and `coi_variant`.

### Running

    coiaf_wrapper(
      snp_calls = "snp_calls.tsv",
      output = "coi_estimates.tsv"
    )

    Rscript exec/coiaf_wrapper \
      --snp_calls snp_calls.tsv \
      --output coi_estimates.tsv

Requires **coiaf** (Suggests).

## See also

[`run_coiaf()`](https://plasmogenepi.github.io/PGEcore/reference/run_coiaf.md),
[`vignette("input-formats", package = "PGEcore")`](https://plasmogenepi.github.io/PGEcore/articles/input-formats.md)
