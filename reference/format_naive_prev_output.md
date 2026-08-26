# Format naive allele-prevalence output for AA or MH input

Format naive allele-prevalence output for AA or MH input

## Usage

``` r
format_naive_prev_output(prevalence, from_aa)
```

## Arguments

- prevalence:

  Prevalence table with `target_name` and `variant`.

- from_aa:

  Logical; `TRUE` when input was amino acid calls.

## Value

Formatted tibble for writing.
