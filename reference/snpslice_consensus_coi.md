# Consensus COI across SNP-Slice restarts

Pools the strain assignments from several independent SNP-Slice restarts
into one per-host complexity of infection (COI) that discounts strains
the restarts do not agree on.

## Usage

``` r
snpslice_consensus_coi(chains, estimate)
```

## Arguments

- chains:

  List of per-restart results, or a single result object.

- estimate:

  Point estimate to read from each restart, `"map"` or `"final_sample"`.

## Value

Numeric vector, one consensus COI per host, floored at 1.

## Why this exists

SNP-Slice fits allele frequencies well but over-parameterizes the strain
dictionary to do so: it adds many low-support strains, often carried by
a single host. Each such strain adds a full +1 to that host's COI while
contributing almost nothing to the frequencies, so the raw row sum of
the allocation matrix over-counts COI. On top of that, restarts are
independent optimizations that land on different dictionaries, so any
single restart's assignments are partly noise. This function addresses
both problems by asking, per host and per strain, how consistently that
assignment is recovered across the different chains, and weighting the
count accordingly.

## How it works

1.  **Match strains across restarts by sequence.** A strain index means
    nothing across restarts, so haplotypes are keyed by their dictionary
    row (the concatenated allele string). Duplicate rows within a
    restart are collapsed before counting.

2.  **Build a membership matrix.** `membership[i, h]` is the fraction of
    restarts in which host `i` carries haplotype `h`. Assignments
    recovered by every restart score 1; those seen in one restart out of
    ten score 0.1.

3.  **Weight each haplotype by cohort support.** Support is
    `colSums(membership)`, the expected number of hosts carrying `h`.
    Each haplotype gets weight `1 - exp(-support / mean(support))`. This
    is a smooth discount rather than a threshold: a strain with no
    support contributes 0, one of average commonness contributes about
    0.63, and a common strain contributes fully. Hard thresholds were
    tested and rejected because the best cutoff differed by population
    and a small wobble in support flipped whole-strain counts.

4.  **Sum and floor.** Host COI is `membership %*% weights`, floored at
    1 so every host is counted as at least one infection.

Normalizing support by its mean, rather than a fixed host count, is what
keeps the estimate from drifting with panel or cohort size. More loci
resolve more strains, which lowers the mean support, so the same
absolute support earns a higher weight. More hosts raise the mean, so
the same support earns a lower weight. Both adjustments are the desired
direction.

In benchmarking on simulated populations this consensus estimate beat
the same weighting applied to a single restart, nearly removed the loci
drift, and was more reproducible between independent seeds. The gain
saturates at roughly three restarts.
