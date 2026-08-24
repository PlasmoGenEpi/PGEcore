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
├── src/               # Install-time C (THEREALMcCOIL)
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
| `estimate_coi_naive` | Pure R | — |
| `estimate_allele_frequency_naive` | Pure R | — |
| `estimate_allele_prevalence_naive` | Pure R | — |
| `allele_per_locus_summary` | Pure R | — |
| `run_coiaf` / `coiaf_wrapper` | R-package wrapper | `coiaf` |
| `filter_biallelic_calls` | Pure R | — |
| `filter_to_highest_diversity_independent_snp_call` | Pure R | — |
| `slaf_from_mhaps_freqs` | Pure R | — |
| `slaf_from_stave_mlaf` | R-package wrapper | `variantstring` |
| `multilocus_prevfreq_naive` | Pure R | — |
| `multilocus_prevfreq_naive_variantstring` | R-package wrapper | `variantstring` |
| `snp_calls_to_vcf` | R-package wrapper | `Biostrings` |
| `vcf_to_snp_calls` | Pure R | — |
| `add_ref_seqs_with_targeted_ref_fasta` | R-package wrapper | `Biostrings` |
| `add_ref_seqs_with_full_genome_ref_fasta` | R-package wrapper | `Biostrings` |
| `pileup_specific_snps` | R-package wrapper | `Biostrings`, `pwalign` |
| `translate_loci_of_interest` | R-package wrapper | `Biostrings`, `pwalign` |
| `per_locus_popgen_summary` | R-package wrapper | `ape`, `msa`, `pegas`; Muscle/Clustal binaries |
| `calculate_fws_from_vcf` | R-package wrapper | `moimix`, `SeqArray` |
| `run_moire` / `moire_wrapper` | R-package wrapper | `moire`, `checkmate` |
| `run_malariaem` / `malariaem_wrapper` | R-package wrapper | `malaria.em`, `checkmate` |
| `dcifer_slaf_wrapper` | R-package wrapper | `dcifer` |
| `dcifer_ibd_wrapper` | R-package wrapper | `dcifer`, `foreach`, `doParallel`, `parallelly`, `iterators` |
| `snpslice_wrapper` | R-package wrapper | `snp.slicer`, `variantstring` |
| `FreqEstimationModel_wrapper` | R-package wrapper | `FreqEstimationModel`, `variantstring`, `foreach`, `doMC`, `plyr`, `coda`, `abind` |
| `IDM_wrapper` | Vendored algorithm wrapper | `Rmpfr`, `openxlsx` |
| `MultiLociBiallelicModel_wrapper` | Vendored algorithm wrapper | `variantstring` |
| `THEREALMcCOIL_wrapper` | Vendored C in `src/` (install-time compile) | — |

## External and optional dependencies

| PGEcore function/script | Dependency | Type | Notes |
| ----------------------- | ---------- | ---- | ----- |
| `count_samples_by_coi` | — | — | Pure R |
| `estimate_coi_naive` | — | — | Pure R |
| `estimate_allele_frequency_naive` | — | — | Pure R |
| `estimate_allele_prevalence_naive` | — | — | Pure R |
| `allele_per_locus_summary` | — | — | Pure R |
| `run_coiaf` / `coiaf_wrapper` | `coiaf` | Optional R (Suggests) | Install separately; not pulled in by PGEcore |
| `filter_biallelic_calls` | — | — | Pure R |
| `filter_to_highest_diversity_independent_snp_call` | — | — | Pure R |
| `slaf_from_mhaps_freqs` | — | — | Pure R |
| `slaf_from_stave_mlaf` | `variantstring` | Optional R (Suggests) | Install separately; not pulled in by PGEcore |
| `multilocus_prevfreq_naive` | — | — | Pure R |
| `multilocus_prevfreq_naive_variantstring` | `variantstring` | Optional R (Suggests) | Install separately; not pulled in by PGEcore |
| `snp_calls_to_vcf` | `Biostrings` | Optional R (Suggests) | Bioconductor; install separately |
| `vcf_to_snp_calls` | — | — | Pure R |
| `add_ref_seqs_with_targeted_ref_fasta` | `Biostrings` | Optional R (Suggests) | Bioconductor; install separately |
| `add_ref_seqs_with_full_genome_ref_fasta` | `Biostrings` | Optional R (Suggests) | Bioconductor; install separately |
| `pileup_specific_snps` | `Biostrings`, `pwalign` | Optional R (Suggests) | Bioconductor; install separately |
| `translate_loci_of_interest` | `Biostrings`, `pwalign` | Optional R (Suggests) | Bioconductor; install separately |
| `per_locus_popgen_summary` | `ape`, `msa`, `pegas` | Optional R (Suggests) | **msa** also needs a `muscle`, `clustalw`, or `clustalo` binary on `PATH` (not installed by the R package) |
| `calculate_fws_from_vcf` | `moimix`, `SeqArray` | Optional R (Suggests) | `SeqArray` is Bioconductor; `moimix` is GitHub-only (`bahlolab/moimix`) |
| `run_moire` / `moire_wrapper` | `moire`, `checkmate` | Optional R (Suggests) | Install separately |
| `run_malariaem` / `malariaem_wrapper` | `malaria.em`, `checkmate` | Optional R (Suggests) | Install separately |
| `dcifer_slaf_wrapper` | `dcifer` | Optional R (Suggests) | Install separately |
| `dcifer_ibd_wrapper` | `dcifer` (+ parallel helpers) | Optional R (Suggests) | `foreach`, `doParallel`, `parallelly`, `iterators` |
| `snpslice_wrapper` | `snp.slicer`, `variantstring` | Optional R (Suggests) | `variantstring` 1.x required |
| `FreqEstimationModel_wrapper` | `FreqEstimationModel` (+ helpers) | Optional R (Suggests) | `variantstring`, `foreach`, `doMC`, `plyr`, `coda`, `abind` |
| `IDM_wrapper` | `Rmpfr`, `openxlsx` | Optional R (Suggests) | Vendored Hashemi & Schneider (2024) Incomplete Data Model |
| `MultiLociBiallelicModel_wrapper` | `variantstring` | Optional R (Suggests) | Vendored SNPModel.R; `variantstring` 1.x required |
| `THEREALMcCOIL_wrapper` | Vendored C (`src/`) | Compiled at install | No runtime `R CMD SHLIB`; grid in `inst/extdata/` |

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
