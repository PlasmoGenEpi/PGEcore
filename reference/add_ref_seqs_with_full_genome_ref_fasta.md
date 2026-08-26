# Add reference sequences extracted from a genome FASTA onto a panel BED table

For each BED interval, extracts `length` bases starting at 0-based
`start` (Biostrings 1-based `start + 1`) from the contig named in
`#chrom`. Minus- strand intervals are reverse-complemented. Genome FASTA
names are truncated at the first whitespace.

## Usage

``` r
add_ref_seqs_with_full_genome_ref_fasta(
  ref_bed,
  genome_fasta,
  output = NULL,
  overwrite = FALSE
)
```

## Arguments

- ref_bed:

  Path to a panel BED TSV, or a data frame with the same columns. See
  *Inputs*.

- genome_fasta:

  Path to a genome FASTA. See *Inputs*.

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

- **`genome_fasta`**: Genome FASTA to extract intervals from.

### Outputs

- **`output`** (optional): `ref_bed` TSV with a `ref_seq` column. If
  `NULL`, results are returned without writing.

### Running

    add_ref_seqs_with_full_genome_ref_fasta(
      ref_bed = "ref_bed.tsv",
      genome_fasta = "genome.fasta",
      output = "ref_bed_with_seq.tsv"
    )

    Rscript exec/add_ref_seqs_with_full_genome_ref_fasta \
      --ref_bed ref_bed.tsv \
      --genome_fasta genome.fasta \
      --output ref_bed_with_seq.tsv

Requires **Biostrings** (Suggests).

## See also

[`vignette("input-formats", package = "PGEcore")`](https://plasmogenepi.github.io/PGEcore/articles/input-formats.md)
