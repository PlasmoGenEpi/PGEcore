# Format SNP calls for the McCOIL categorical model

Format SNP calls for the McCOIL categorical model

## Usage

``` r
prep_input_categorical(df)
```

## Arguments

- df:

  Output of
  [`read_and_preprocess_snp_call()`](https://plasmogenepi.github.io/PGEcore/reference/read_and_preprocess_snp_call.md).

## Value

A data frame (samples x sites) of scores `1` / `0` / `0.5` / `-1`.
