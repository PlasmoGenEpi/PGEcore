# Getting started with PGEcore

PGEcore is an R package for malaria genomics analysis. It **wraps
existing tools** so they share a consistent R API, CLI, and TSV formats,
and it also **adds extra pieces**—for example naive COI, frequency, and
prevalence methods, filters, and converters—that sit alongside those
wrappers.

Every tool is available both as an **R function** and as a
**command-line program**, so you can use it interactively, in scripts,
or as a standard library across workflow pipelines. Shared table layouts
make it easier to move results from one analysis into the next without
reformatting.

## Install

``` r

# install.packages("pak")  # if you do not have it yet
# pak::pak("PlasmoGenEpi/PGEcore")
library(PGEcore)
```

Specialised packages used by some wrappers (for example `dcifer`,
`moire`) are **Suggests** and must be installed separately when you need
those tools.

## From R

``` r

library(PGEcore)

# File paths
count_samples_by_coi("my_coi_table.tsv", output = "coi_distribution.tsv")

# Or in-memory tables
coi_table <- readr::read_tsv("my_coi_table.tsv", show_col_types = FALSE)
count_samples_by_coi(coi_table)
```

Function names match the CLI name (`count_samples_by_coi`,
`moire_wrapper`, …). That is the API to start with. Each tool’s help
page has the same structure—**Inputs**, **Outputs**, and **Running**
(R + CLI examples). Shared table layouts live in
[`vignette("input-formats", package = "PGEcore")`](https://plasmogenepi.github.io/PGEcore/articles/input-formats.md).

A small number of wrappers also export an in-memory `run_*` helper
(`run_coiaf`, `run_moire`, `run_malariaem`) for when you already have
data frames or model objects in R; see
[`?run_moire`](https://plasmogenepi.github.io/PGEcore/reference/run_moire.md).

``` r

?count_samples_by_coi
help(package = "PGEcore")
```

## Command line

Each tool has an executable in the package `exec/` directory. Installing
the package does not put those on `PATH`, so ask R where they landed and
add that directory yourself:

``` bash
PGECORE_EXEC="$(Rscript -e 'cat(system.file("exec", package="PGEcore"))')"
export PATH="$PGECORE_EXEC:$PATH"

count_samples_by_coi \
  --coi_table my_coi_table.tsv \
  --output coi_distribution.tsv
```

That `export` lasts only for the current shell; put both lines in your
`~/.bashrc` or `~/.zshrc` to make it permanent. Capturing the location
in a variable first is deliberate:
[`system.file()`](https://rdrr.io/r/base/system.file.html) returns an
empty string when the package cannot be found, and an empty entry in
`PATH` means the current directory. Check that `PGECORE_EXEC` actually
holds a path; if it is empty, the package most likely isn’t installed
yet.

Or use the full path directly, without touching `PATH`:

``` bash
"$PGECORE_EXEC/count_samples_by_coi" \
  --coi_table my_coi_table.tsv \
  --output coi_distribution.tsv
```

From a source checkout of this repository.
[`pak::local_install()`](https://pak.r-lib.org/reference/local_install.html)
also installs the package’s dependencies, which a plain
`R CMD INSTALL .` does not:

``` bash
git clone https://github.com/PlasmoGenEpi/PGEcore.git
cd PGEcore
Rscript -e 'pak::local_install(".")'

exec/count_samples_by_coi \
  --coi_table inst/extdata/example_coi_table.tsv \
  --output coi_distribution.tsv

exec/estimate_coi_naive \
  --allele_table inst/extdata/example_allele_table.tsv \
  --output coi_table.tsv

exec/dcifer_slaf_wrapper \
  --allele_table inst/extdata/example_allele_table.tsv \
  --slaf_output slaf.tsv
```

Pass ordinary file paths. Example TSVs that show the expected columns
are under `inst/extdata/` (see the **input-formats** vignette).

## Combining analyses

1.  Match your tables to a **standard format** (SNP calls, allele table,
    amino acid calls, COI table, …).
2.  Run a PGEcore function or CLI.
3.  Feed the output TSV into the next tool when the column layouts
    already align.

That interoperability—one invocation style and shared table layouts—is a
core goal of the package.
