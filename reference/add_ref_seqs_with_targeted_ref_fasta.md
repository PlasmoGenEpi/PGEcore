# Add reference sequences from a targeted FASTA onto a panel BED table

Joins `ref_seq` onto a panel location table by matching FASTA record
names to `target_name`.

## Usage

``` r
add_ref_seqs_with_targeted_ref_fasta(
  ref_bed,
  target_fasta,
  output = NULL,
  overwrite = FALSE
)
```

## Arguments

- ref_bed:

  Path to a panel BED TSV, or a data frame with the same columns. See
  *Inputs*.

- target_fasta:

  Path to a FASTA whose record names match `target_name`.

- output:

  Optional output TSV path.

- overwrite:

  If `FALSE` (default), refuse to overwrite `output`.

## Value

A tibble of `ref_bed` with a `ref_seq` column.

## Details

### Inputs

- **`ref_bed`**: Panel BED TSV with header (`#chrom`, `start`, `end`,
  `target_name`, `length`, `strand`), as a file path or data frame.

- **`target_fasta`**: FASTA whose record names match `target_name`.

### Outputs

- **`output`** (optional): `ref_bed` TSV with a `ref_seq` column. If
  `NULL`, results are returned without writing.

### Running

    add_ref_seqs_with_targeted_ref_fasta(
      ref_bed = "ref_bed.tsv",
      target_fasta = "targets.fasta",
      output = "ref_bed_with_seq.tsv"
    )

    Rscript exec/add_ref_seqs_with_targeted_ref_fasta \
      --ref_bed ref_bed.tsv \
      --target_fasta targets.fasta \
      --output ref_bed_with_seq.tsv

Requires **Biostrings** (Suggests) to read the FASTA.

## See also

[`vignette("input-formats", package = "PGEcore")`](https://plasmogenepi.github.io/PGEcore/articles/input-formats.md)
