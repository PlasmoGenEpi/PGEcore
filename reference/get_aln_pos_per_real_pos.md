# Map a gapless sequence coordinate to an aligned coordinate

Map a gapless sequence coordinate to an aligned coordinate

## Usage

``` r
get_aln_pos_per_real_pos(aligned_seq, pos)
```

## Arguments

- aligned_seq:

  Aligned sequence that may contain `"-"` characters.

- pos:

  1-based position in the gapless sequence.

## Value

The corresponding 1-based index in `aligned_seq`.
