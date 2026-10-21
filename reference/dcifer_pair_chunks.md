# Split pair indices into one non-empty chunk per worker

Split pair indices into one non-empty chunk per worker

## Usage

``` r
dcifer_pair_chunks(total_pairs, total_cores)
```

## Arguments

- total_pairs:

  Number of sample pairs.

- total_cores:

  Number of workers.

## Value

A list of integer index vectors covering `seq_len(total_pairs)` exactly
once.
