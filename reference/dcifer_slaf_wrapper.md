# Estimate single-locus allele frequencies with Dcifer

Estimates per-locus allele frequencies from an allele table (optionally
with a COI table). Requires **dcifer** (Suggests).

## Usage

``` r
dcifer_slaf_wrapper(
  allele_table,
  slaf_output,
  coi_table = NULL,
  specimen_name_col = "specimen_name",
  target_name_col = "target_name",
  target_value_col = "seq",
  tol = 1e-04,
  qstart = 0.5,
  coi_lrank = 2L
)
```

## Arguments

- allele_table:

  Path to allele table TSV. See *Inputs*.

- slaf_output:

  Path for SLAF TSV output. See *Outputs*.

- coi_table:

  Optional path to COI table TSV. See *Inputs*.

- specimen_name_col, target_name_col, target_value_col:

  Column names in `allele_table` / `coi_table`.

- tol, qstart:

  Passed to
  [`dcifer::calcAfreq()`](https://eppicenter.github.io/dcifer/reference/calcAfreq.html).

- coi_lrank:

  Rank of the locus used by
  [`dcifer::getCOI()`](https://eppicenter.github.io/dcifer/reference/getCOI.html)
  when `coi_table` is not supplied.

## Value

The SLAF tibble (also written to `slaf_output`).

## Details

### Inputs

- **`allele_table`**: Allele table TSV. See
  [`vignette("input-formats", package = "PGEcore")`](https://plasmogenepi.github.io/PGEcore/articles/input-formats.md).

- **`coi_table`**: Optional COI table TSV. If omitted, COI is inferred
  with
  [`dcifer::getCOI()`](https://eppicenter.github.io/dcifer/reference/getCOI.html).

### Outputs

- **`slaf_output`**: Single-locus allele frequencies (default columns
  `target_name`, `seq`, `freq`, `sample_total`).

### Running

    dcifer_slaf_wrapper(
      allele_table = "allele_table.tsv",
      slaf_output = "slaf.tsv"
    )

    Rscript exec/dcifer_slaf_wrapper \
      --allele_table allele_table.tsv \
      --slaf_output slaf.tsv

Requires **dcifer** (Suggests).

## See also

[`vignette("input-formats", package = "PGEcore")`](https://plasmogenepi.github.io/PGEcore/articles/input-formats.md),
[`dcifer_ibd_wrapper()`](https://plasmogenepi.github.io/PGEcore/reference/dcifer_ibd_wrapper.md)
