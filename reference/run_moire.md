# Run MOIRE MCMC analysis

Runs MOIRE MCMC on a prepared `moire_object`. Requires **moire**
(Suggests). For reading allele tables and writing summary TSVs, use
[`moire_wrapper()`](https://plasmogenepi.github.io/PGEcore/reference/moire_wrapper.md).

## Usage

``` r
run_moire(moire_object)
```

## Arguments

- moire_object:

  List with `moire_data` and `moire_parameters`. See *Inputs*.

## Value

The object returned by
[`moire::run_mcmc()`](https://EPPIcenter.github.io/moire/reference/run_mcmc.html).

## Details

### Inputs

- **`moire_object`**: List with `moire_data` and `moire_parameters`, as
  created by
  [`create_moire_input()`](https://plasmogenepi.github.io/PGEcore/reference/create_moire_input.md)
  (via
  [`moire_wrapper()`](https://plasmogenepi.github.io/PGEcore/reference/moire_wrapper.md)).

### Outputs

- Returns the object from
  [`moire::run_mcmc()`](https://EPPIcenter.github.io/moire/reference/run_mcmc.html)
  (not written to disk).

### Running

    run_moire(moire_object)

File and CLI users should call
[`moire_wrapper()`](https://plasmogenepi.github.io/PGEcore/reference/moire_wrapper.md)
/ `Rscript exec/moire_wrapper ...`.

Requires **moire** (Suggests).

## See also

[`moire_wrapper()`](https://plasmogenepi.github.io/PGEcore/reference/moire_wrapper.md),
[`vignette("input-formats", package = "PGEcore")`](https://plasmogenepi.github.io/PGEcore/articles/input-formats.md)
