# PGEcore

An R package for malaria genomics analysis. PGEcore does two things:

1. **Wraps existing tools** (for example `coiaf`, `moire`, `dcifer`) behind a
   consistent R API and CLI, with shared TSV inputs and outputs so steps are
   easy to chain.
2. **Adds extra analyses** that are not just thin wrappers—naive COI /
   frequency / prevalence estimators, filters, format converters, summaries,
   and similar helpers.

Analyses are available both as **R functions** and as **command-line tools**,
so you can use PGEcore interactively, in your own scripts, or as a shared
library across workflow pipelines (for example Nextflow or WDL).

> **PGEcore does not install the specialised software it wraps.** Optional
> dependencies (for example the `coiaf` package) must be installed separately
> when you need those tools.

## Install

```r
# From GitHub
# remotes::install_github("PlasmoGenEpi/PGEcore")

# From a local clone
devtools::install()
```

## Two ways to run every tool

1. An **R function** — interactive analysis, notebooks, or other R code
2. A **CLI** under `exec/` — shell scripts or batch jobs

Flags and file formats match between the two.

### From R

```r
library(PGEcore)

count_samples_by_coi("my_coi_calls.tsv", output = "coi_distribution.tsv")

# Or pass data frames already in memory
# run_coiaf(snp_data = snp_df)   # requires the coiaf package
```

```r
help(package = "PGEcore")
browseVignettes("PGEcore")
# vignette("getting-started", package = "PGEcore")
# vignette("input-formats", package = "PGEcore")
```

### Command line

Once the package executables are on your `PATH` (for example after a Conda
install of `r-pgecore`):

```bash
count_samples_by_coi \
  --coi_calls my_coi_calls.tsv \
  --output coi_distribution.tsv
```

From a source checkout of this repository:

```bash
Rscript exec/count_samples_by_coi \
  --coi_calls inst/extdata/example_coi_table.tsv \
  --output coi_distribution.tsv

Rscript exec/coiaf_wrapper \
  --snp_data inst/extdata/example_collapsed_snp_calls.tsv \
  --output coi_estimates.tsv
```

Pass ordinary file paths to your data. Bundled examples for trying formats live
in `inst/extdata/`.

## Standard input formats

Tools share a small set of TSV layouts so outputs from one step can feed the
next. Required columns (minimum):

| Format | Required columns | Example file |
| ------ | ---------------- | ------------ |
| COI calls | `specimen_name`, `coi` | `inst/extdata/example_coi_table.tsv` |
| SNP calls | `specimen_name`, `snp_name`, `reads`, `seq_base` | `inst/extdata/example_collapsed_snp_calls.tsv` |
| Allele / microhap table | `specimen_name`, `target_name`, `seq`, `reads` | `inst/extdata/example_allele_table.tsv` |
| Amino acid calls | `specimen_name`, `gene_id`, `aa_position`, `aa`, `reads` (+ often `target_name`, `aa_locus`, …) | `inst/extdata/example_amino_acid_calls.tsv` |
| Loci groups | `group_id`, `gene_id`, `aa_position` | `inst/extdata/example_loci_groups.tsv` |
| Allele frequency (MLAF) | `group_id`, `variant`, `freq` | `inst/extdata/example_mlaf.tsv` |
| Population MAF (PLMAF) | `snp_name`, `seq_base`, `plmaf` | `inst/extdata/example_coiaf_plmaf.tsv` |

Full column notes and which tools consume each format: vignette **`input-formats`**.

## Available tools

| CLI / function | Kind | Optional dependency |
| -------------- | ---- | ------------------- |
| `count_samples_by_coi` | Pure R | — |
| `estimate_coi_naive` | Pure R | — |
| `estimate_allele_frequency_naive` | Pure R | — |
| `estimate_allele_prevalence_naive` | Pure R | — |
| `allele_per_locus_summary` | Pure R | — |
| `coiaf_wrapper` / `run_coiaf` | Wrapper | `coiaf` |
| `filter_biallelic_calls` | Pure R | — |
| `filter_to_highest_diversity_independent_snp_call` | Pure R | — |
| `slaf_from_mhaps_freqs` | Pure R | — |
| `slaf_from_stave_mlaf` | Wrapper | `variantstring` |
| `multilocus_prevfreq_naive` | Pure R | — |
| `multilocus_prevfreq_naive_variantstring` | Wrapper | `variantstring` |
| `snp_calls_to_vcf` | Wrapper | `Biostrings` |
| `vcf_to_snp_calls` | Pure R | — |
| `add_ref_seqs_with_targeted_ref_fasta` | Wrapper | `Biostrings` |
| `add_ref_seqs_with_full_genome_ref_fasta` | Wrapper | `Biostrings` |
| `pileup_specific_snps` | Wrapper | `Biostrings`, `pwalign` |
| `translate_loci_of_interest` | Wrapper | `Biostrings`, `pwalign` |
| `per_locus_popgen_summary` | Wrapper | `ape`, `msa`, `pegas` (+ Muscle/Clustal on PATH) |
| `calculate_fws_from_vcf` | Wrapper | `moimix`, `SeqArray` |
| `moire_wrapper` / `run_moire` | Wrapper | `moire` |
| `malariaem_wrapper` / `run_malariaem` | Wrapper | `malaria.em` |
| `dcifer_slaf_wrapper` | Wrapper | `dcifer` |
| `dcifer_ibd_wrapper` | Wrapper | `dcifer` (+ parallel helpers) |
| `snpslice_wrapper` | Wrapper | `snp.slicer`, `variantstring` |
| `FreqEstimationModel_wrapper` | Wrapper | `FreqEstimationModel` (+ helpers) |
| `IDM_wrapper` | Vendored algorithm | `Rmpfr`, `openxlsx` |
| `MultiLociBiallelicModel_wrapper` | Vendored algorithm | `variantstring` |
| `THEREALMcCOIL_wrapper` | Vendored C (`src/`) | — |

## Optional dependencies

Wrappers that call another R package list that package under **Suggests**.
Install only what you need, for example:

```r
install.packages(
  "coiaf",
  repos = c("https://plasmogenepi.r-universe.dev", "https://cloud.r-project.org")
)
```

With Conda, a minimal environment might look like:

```yaml
dependencies:
  - r-base
  - r-pgecore
  - r-coiaf    # only if you use coiaf_wrapper
```

## Package layout

```text
PGEcore/
├── R/              # Exported API and helpers
├── exec/           # Thin CLIs (same flags as the R API)
├── src/            # THEREALMcCOIL C (compiled at install)
├── inst/extdata/   # Example inputs (standard formats)
├── vignettes/      # Getting started + input formats
├── man/
└── tests/
```

## Contribute

PRs go to `develop` (Gitflow). To add a tool:

1. Implement in `R/` (validate → prepare → run → format).
2. Export one high-level function; keep helpers internal.
3. Put specialised deps in `Suggests` and use `check_suggested_pkg()`.
4. Add a thin `exec/` CLI with the same optparse flags.
5. Add tests and document formats in the `input-formats` vignette if you
   introduce a new shared table layout.
