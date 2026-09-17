# PGEcore 0.1.0

* R package for malaria genomics analysis (`R/`, `exec/`, `tests/`,
  `inst/extdata/`, `vignettes/`).
* Dual interface: R API and command-line tools with matching flags and shared
  TSV formats.
* Vignettes: `getting-started` and `input-formats`.
* Specialised analysis packages remain **Suggests** — installing PGEcore does
  not pull them in.
* `THEREALMcCOIL` C is compiled at install time from `src/`.
