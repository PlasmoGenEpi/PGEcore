# PGEcore 0.1.0

* R package for malaria genomics analysis (`R/`, `exec/`, `tests/`,
  `inst/extdata/`, `vignettes/`).
* Dual interface: R API and command-line tools with matching flags and shared
  TSV formats.
* Vignettes: `getting-started` and `input-formats`.
* Specialised analysis packages remain **Suggests** — installing PGEcore does
  not pull them in.
* `THEREALMcCOIL` C is compiled at install time from `src/`.
* `THEREALMcCOIL_wrapper(model = "proportional")` gives each sample one set of
  strain proportions shared by all its loci, estimates read noise, and treats
  minor alleles below the detection limit as censored. The `epsilon`
  argument and the fitted Beta grid it used are removed. COI and
  allele-frequency outputs are medians pooled over all chains, and `threads`
  sets how many chains run at once. `M0` takes one starting COI per chain
  (default `"1,5,15"`) so R-hat can reveal a chain that has not converged.
* McCOIL error-rate proposals used a variance as the standard deviation; fixed.
  `err_method = 2` now refreshes its likelihoods after each redraw. The wrapper
  rejects `err_method = 2`, and any `err_method` other than 1 for the
  proportional model.
* `filter_to_highest_diversity_independent_snp_call()` gains
  `drop_artifact_snps` (drop SNPs whose minor allele recurs in otherwise
  monoclonal specimens, such as PCR stutter) and `exclude_snp_names`; both act
  before SNPs are ranked.
