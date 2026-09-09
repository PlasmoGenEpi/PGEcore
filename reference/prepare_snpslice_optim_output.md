# Format SNP-Slice per-restart optimization diagnostics

One row per restart: `chain_id`, `seed`, `map_logpost`, `is_best`,
`map_iteration`, `final_iteration`, `plateau_frac` (`map_iteration`
divided by `final_iteration`), `map_kstar`, `map_ktrunc`, `coi_mean`,
and `coi_ccc_to_best` (Lin's CCC between that restart's per-host COI and
the reported restart's).

## Usage

``` r
prepare_snpslice_optim_output(snpslice_res)
```
