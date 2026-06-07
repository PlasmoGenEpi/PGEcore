#!/usr/bin/env Rscript
# calculate_fws_from_vcf.R
# -----------------------------------------------------------------------------
# Convert a VCF to GDS (only when needed) and calculate the Fws within-host
# diversity statistic for each sample using moimix::getFws().
#
# Example:
#   ./calculate_fws_from_vcf.R -i input.vcf.gz -o fws_result.tsv
# -----------------------------------------------------------------------------
#setRepositories(ind = 1:3);
#install.packages(c("devtools", "remotes", "SeqArray", "optparse", "dplyr", "readr", "tidyr"))
#install.packages(c('SeqVarTools', 'BiocParallel'))
#devtools::install_github("bahlolab/moimix")
suppressPackageStartupMessages({
  library(optparse)
  library(moimix)
  library(SeqArray)
  library(dplyr)
  library(tidyr)
  library(readr)
})
# ---- Argument parsing --------------------------------------------------------
option_list <- list(
  make_option(c("-i", "--input"), type = "character", default = NULL,
              metavar = "FILE",
              help = "Input VCF file (.vcf or .vcf.gz) [required]"),
  make_option(c("-o", "--output"), type = "character", default = "fws_result.tsv",
              metavar = "FILE",
              help = "Output TSV of Fws results [default: %default]"),
  make_option(c("-g", "--gds"), type = "character", default = NULL,
              metavar = "FILE",
              help = paste("GDS file path. If omitted, it is derived from the",
                           "input VCF by replacing the .vcf/.vcf.gz suffix with .gds")),
  make_option(c("-p", "--population_name"), type = "character", default = NULL,
              metavar = "NAME",
              help = paste("Optional population name. If provided, adds a",
                           "'population_name' column with this value to the output",
                           "(useful for combining multiple populations later)")),
  make_option(c("-f", "--force"), action = "store_true", default = FALSE,
              help = "Force re-creation of the GDS even if an up-to-date one exists"),
  make_option(c("-v", "--verbose"), action = "store_true", default = FALSE,
              help = "Print progress messages (otherwise run silently)")
)
parser <- OptionParser(
  usage = "%prog -i input.vcf.gz -o fws_result.tsv [options]",
  option_list = option_list,
  description = "Calculate Fws from a VCF using moimix::getFws()."
)
opt <- parse_args(parser)
# Emit a message only when --verbose is set.
say <- function(...) { if (isTRUE(opt$verbose)) message(sprintf(...)) }
# ---- Validate input ----------------------------------------------------------
if (is.null(opt$input)) {
  print_help(parser)
  stop("An input VCF (-i/--input) is required.", call. = FALSE)
}
if (!file.exists(opt$input)) {
  stop(sprintf("Input VCF not found: %s", opt$input), call. = FALSE)
}
# Derive the GDS path from the VCF name if one wasn't supplied.
if (is.null(opt$gds)) {
  opt$gds <- sub("\\.vcf(\\.gz)?$", ".gds", opt$input, ignore.case = TRUE)
  if (identical(opt$gds, opt$input)) {
    # Input didn't end in .vcf/.vcf.gz; just append .gds.
    opt$gds <- paste0(opt$input, ".gds")
  }
}
# ---- Decide whether the GDS needs to be (re)created --------------------------
needs_conversion <- TRUE
if (!opt$force && file.exists(opt$gds)) {
  if (file.mtime(opt$gds) >= file.mtime(opt$input)) {
    needs_conversion <- FALSE
    say("Using existing GDS (up to date with VCF): %s", opt$gds)
  } else {
    say("Existing GDS is older than the VCF; re-creating: %s", opt$gds)
  }
} else if (opt$force && file.exists(opt$gds)) {
  say("--force set; re-creating existing GDS: %s", opt$gds)
}
if (needs_conversion) {
  say("Converting VCF -> GDS: %s -> %s", opt$input, opt$gds)
  if (isTRUE(opt$verbose)) {
    seqVCF2GDS(opt$input, opt$gds)
  } else {
    suppressMessages(seqVCF2GDS(opt$input, opt$gds, verbose = FALSE))
  }
}
# ---- Calculate Fws -----------------------------------------------------------
gds <- seqOpen(opt$gds)
fws_result <- if (isTRUE(opt$verbose)) getFws(gds) else suppressMessages(getFws(gds))
seqClose(gds)
fws_result_df <- tibble(
  specimen_name = names(fws_result),
  fws = as.numeric(fws_result)
) %>%
  arrange(fws)
# Optionally tag every row with the population name for later merging.
if (!is.null(opt$population_name)) {
  fws_result_df <- fws_result_df %>%
    mutate(population_name = opt$population_name)
}
write_tsv(fws_result_df, opt$output)
say("Wrote %d samples to %s", nrow(fws_result_df), opt$output)
