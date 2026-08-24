# PGEcore

Shared R functions and command-line interfaces used across PlasmoGenEpi
Nextflow and WDL pipelines.

> **PGEcore does not install or distribute the external software or specialised
> R packages that it wraps.** Those tools must be installed separately and
> available at runtime (for example on `PATH`, or as Conda packages in the same
> process environment). This allows each workflow process to install only the
> software it needs.

This repository is being converted from a Git-submodule script collection into a
versioned R package. Legacy modules remain under `scripts/` while they are
migrated into `R/` + `exec/`.

## Install

```r
# From GitHub (development)
# remotes::install_github("PlasmoGenEpi/PGEcore")

# Or from a local clone
devtools::install()
```

Specialised Suggests packages (for example `coiaf`) are **not** installed
automatically. Install them separately when you need a given wrapper.

## Two ways to use PGEcore

### From R

```r
library(PGEcore)

coi_path <- system.file("extdata", "example_coi_table.tsv", package = "PGEcore")
count_samples_by_coi(coi_path)

# Optional dependency (install coiaf separately)
# results <- run_coiaf(snp_data = snp_df, plmaf = NULL)
```

### From a workflow / CLI

After installation, thin CLIs live in the package `exec/` directory (executable
scripts shipped with the package). Conda recipes for `r-pgecore` should expose
them on `PATH`; until then, invoke via the installed path:

```bash
# Resolve installed CLI path
COI_CLI="$(Rscript -e 'cat(system.file("exec", "count_samples_by_coi", package = "PGEcore"))')"
"$COI_CLI" --coi_calls example_coi_table.tsv --output coi_distribution.tsv

# Or, once Conda puts them on PATH:
count_samples_by_coi --coi_calls example_coi_table.tsv --output coi_distribution.tsv

# Requires the coiaf R package in the same environment
COIAF_CLI="$(Rscript -e 'cat(system.file("exec", "coiaf_wrapper", package = "PGEcore"))')"
"$COIAF_CLI" --snp_data snps.tsv --output coi_estimates.tsv
```

Preserve existing flag names when migrating pipelines from submodule
`Rscript scripts/...` invocations; only the executable path should change.

## Package layout

```text
PGEcore/
├── DESCRIPTION
├── R/                 # Implementation and exported API
├── exec/              # Thin CLIs installed onto PATH
├── inst/extdata/      # Example inputs
├── man/
├── tests/testthat/
├── scripts/           # Legacy modules (migration in progress)
└── data/              # Legacy examples (also in inst/extdata)
```

## Migrated so far

| Function / CLI | Kind | Optional dependency |
| -------------- | ---- | ------------------- |
| `count_samples_by_coi` | Pure R | — |
| `run_coiaf` / `coiaf_wrapper` | R-package wrapper | `coiaf` |

## External and optional dependencies

| PGEcore function/script | Dependency | Type | Notes |
| ----------------------- | ---------- | ---- | ----- |
| `count_samples_by_coi` | — | — | Pure R |
| `run_coiaf` / `coiaf_wrapper` | `coiaf` | Optional R (Suggests) | Install separately; not pulled in by PGEcore |

R packages used across most of PGEcore (see `Imports` in `DESCRIPTION`) are
installed with the package. Specialised analysis packages belong in `Suggests`
and are checked at runtime with `check_suggested_pkg()`.

## Conda / workflow note

A process environment should look conceptually like:

```yaml
dependencies:
  - r-base
  - r-pgecore
  - r-coiaf          # only if this process needs coiaf_wrapper
```

Do **not** build a single environment containing every tool PGEcore can wrap.

## How to contribute

We use [Gitflow](https://www.atlassian.com/git/tutorials/comparing-workflows/gitflow-workflow).
Open PRs into `develop` (not `main`).

### Adding a new wrapper (package style)

1. Put implementation functions in `R/<name>.R` (validate → prepare → run → format).
2. Export only the high-level API; keep helpers internal.
3. If the wrapper needs a specialised R package, list it under `Suggests` and call
   `check_suggested_pkg()` before use. Do not add it to `Imports`.
4. If it shells out to a binary, call `check_external_tool()` and document the
   executable; do not download or install it from PGEcore.
5. Add a thin CLI under `exec/` that only parses arguments (optparse) and calls
   the exported function. Preserve existing flag names when migrating a legacy
   script.
6. Add unit tests that do not require optional tools; integration tests should
   `skip_if_not_installed()` / skip when binaries are missing.
7. Document with roxygen2 and update this README’s dependency table.

Legacy contribution notes for unmigrated `scripts/` modules still apply until
those modules are moved; see per-script READMEs under `scripts/`.
