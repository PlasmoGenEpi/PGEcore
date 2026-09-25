# PGEcore 0.1.0

* R package for malaria genomics analysis (`R/`, `exec/`, `tests/`,
  `inst/extdata/`, `vignettes/`).
* Dual interface: R API and command-line tools with matching flags and shared
  TSV formats.
* Vignettes: `getting-started` and `input-formats`.
* Specialised analysis packages remain **Suggests** — installing PGEcore does
  not pull them in.
* `THEREALMcCOIL` C is compiled at install time from `src/`.
* `THEREALMcCOIL_wrapper()` gains `model = "proportional_joint"`: the
  proportional model with each sample's strain proportions shared across loci
  (the original model integrates over them separately at every locus), read
  noise estimated, and minor alleles below the detection limit treated as
  censored; its chains start spread from 1 to 15 strains so R-hat can reveal a
  stuck chain. COI and allele-frequency outputs are now medians pooled over
  all chains, and `threads` sets how many chains run at once.
* McCOIL error-rate proposals used a variance as the standard deviation; fixed.
  `err_method = 2` now refreshes its likelihoods after each redraw. The wrapper
  rejects `err_method = 2`, and `err_method = 3` for the proportional models.
* `filter_to_highest_diversity_independent_snp_call()` gains
  `drop_artifact_snps` (drop SNPs whose minor allele recurs in otherwise
  monoclonal specimens, such as PCR stutter) and `exclude_snp_names`; both act
  before SNPs are ranked.
