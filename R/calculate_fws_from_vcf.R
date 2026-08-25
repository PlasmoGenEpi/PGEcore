#' Derive a GDS path from a VCF path
#'
#' @param vcf_path Input VCF path.
#' @param gds Optional explicit GDS path.
#' @return GDS path to use.
#' @keywords internal
derive_gds_path <- function(vcf_path, gds = NULL) {
  if (!is.null(gds) && nzchar(gds)) {
    return(gds)
  }
  derived <- sub("\\.vcf(\\.gz)?$", ".gds", vcf_path, ignore.case = TRUE)
  if (identical(derived, vcf_path)) {
    paste0(vcf_path, ".gds")
  } else {
    derived
  }
}

#' Whether a VCF should be (re)converted to GDS
#'
#' @param vcf_path Input VCF path.
#' @param gds_path GDS path.
#' @param overwrite If `TRUE`, always convert.
#' @return Logical.
#' @keywords internal
gds_needs_conversion <- function(vcf_path, gds_path, overwrite = FALSE) {
  if (isTRUE(overwrite)) {
    return(TRUE)
  }
  if (!file.exists(gds_path)) {
    return(TRUE)
  }
  file.mtime(gds_path) < file.mtime(vcf_path)
}

#' Calculate within-host Fws from a VCF via moimix
#'
#' Converts the VCF to GDS with **SeqArray** when the GDS is missing, older
#' than the VCF, or `overwrite` is `TRUE`, then runs `moimix::getFws()`.
#' Requires **moimix** and **SeqArray** (Suggests). `moimix` is installed from
#' GitHub (`bahlolab/moimix`), not CRAN. The VCF must carry per-sample allelic
#' depths (`FORMAT/AD`).
#'
#' @param vcf Input VCF path (`.vcf` or `.vcf.gz`).
#' @param output Output TSV path. Defaults to `"fws_result.tsv"`.
#' @param gds Optional GDS path. If `NULL`, derived by replacing `.vcf` /
#'   `.vcf.gz` with `.gds`.
#' @param population_name Optional population label added as a column.
#' @param overwrite If `TRUE`, rebuild the GDS even when it is up to date.
#' @param verbose If `TRUE`, print progress messages.
#'
#' @return A tibble with `specimen_name`, `fws`, and optionally
#'   `population_name`, sorted by `fws`.
#'
#' @export
calculate_fws_from_vcf <- function(vcf,
                                   output = "fws_result.tsv",
                                   gds = NULL,
                                   population_name = NULL,
                                   overwrite = FALSE,
                                   verbose = FALSE) {
  check_suggested_pkg("SeqArray", "calculate_fws_from_vcf()")
  check_suggested_pkg("moimix", "calculate_fws_from_vcf()")
  if (is.null(vcf) || !nzchar(vcf)) {
    stop("An input VCF (-i/--vcf) is required.", call. = FALSE)
  }
  if (!file.exists(vcf)) {
    stop(sprintf("Input VCF not found: %s", vcf), call. = FALSE)
  }

  say <- function(...) {
    if (isTRUE(verbose)) {
      message(sprintf(...))
    }
  }

  gds_path <- derive_gds_path(vcf, gds)
  needs_conversion <- gds_needs_conversion(vcf, gds_path, overwrite)
  if (!needs_conversion) {
    say("Using existing GDS (up to date with VCF): %s", gds_path)
  } else if (isTRUE(overwrite) && file.exists(gds_path)) {
    say("--overwrite set; re-creating existing GDS: %s", gds_path)
  } else if (file.exists(gds_path)) {
    say("Existing GDS is older than the VCF; re-creating: %s", gds_path)
  }

  if (needs_conversion) {
    say("Converting VCF -> GDS: %s -> %s", vcf, gds_path)
    if (isTRUE(verbose)) {
      SeqArray::seqVCF2GDS(vcf, gds_path)
    } else {
      suppressMessages(SeqArray::seqVCF2GDS(vcf, gds_path, verbose = FALSE))
    }
  }

  gds_obj <- SeqArray::seqOpen(gds_path)
  on.exit(SeqArray::seqClose(gds_obj), add = TRUE)
  fws_result <- if (isTRUE(verbose)) {
    moimix::getFws(gds_obj)
  } else {
    suppressMessages(moimix::getFws(gds_obj))
  }

  fws_result_df <- tibble::tibble(
    specimen_name = names(fws_result),
    fws = as.numeric(fws_result)
  ) |>
    dplyr::arrange(.data$fws)

  if (!is.null(population_name) && nzchar(population_name)) {
    fws_result_df <- fws_result_df |>
      dplyr::mutate(population_name = population_name)
  }
  readr::write_tsv(fws_result_df, output)
  say("Wrote %d samples to %s", nrow(fws_result_df), output)
  fws_result_df
}
