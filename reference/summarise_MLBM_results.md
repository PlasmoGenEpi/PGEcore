# Summarise MLBM haplotype frequencies as variant strings

Summarise MLBM haplotype frequencies as variant strings

## Usage

``` r
summarise_MLBM_results(MLBM_res, MLBM_object, group_name)
```

## Arguments

- MLBM_res:

  Result of
  [`run_MultiLociBiallelicModel()`](https://plasmogenepi.github.io/PGEcore/reference/run_MultiLociBiallelicModel.md).

- MLBM_object:

  Result of
  [`create_MultiLociBiallelicModel_input()`](https://plasmogenepi.github.io/PGEcore/reference/create_MultiLociBiallelicModel_input.md).

- group_name:

  Group id matching `MLBM_object$by_group_table`.

## Value

Tibble with `group_id`, `variant`, `freq`.
