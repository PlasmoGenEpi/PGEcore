# PGEcore 0.1.0

* Converted the repository into an R package (`DESCRIPTION`, `R/`, `exec/`,
  `tests/`, `inst/extdata/`).
* Shared helpers: `check_suggested_pkg()`, `check_external_tool()`, CLI/IO
  validation utilities, and `convert_single_locus_table_to_stave()`.
* Migrated essentially all legacy `scripts/` modules to `R/` + thin `exec/`
  CLIs, preserving optparse flag names where practical.
* `THEREALMcCOIL` C (`McCOIL_categorical`, `McCOIL_prop`) is compiled at
  install time from `src/` and registered for `.C(..., PACKAGE = "PGEcore")`;
  the proportional-method beta grid ships in `inst/extdata/`.
* Specialised analysis packages and PATH tools remain **Suggests** / runtime
  checks — installing PGEcore does not pull them in.

## Known deferrals

* Legacy `scripts/` copies are retained during the transition so submodule
  pipelines keep working; prefer the package API / `exec/` CLIs going forward.
