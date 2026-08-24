#' Add reference sequences from a targeted FASTA onto a panel BED table
#'
#' Joins `ref_seq` onto a panel location table by matching FASTA record names
#' to `target_name`. Requires **Biostrings** (Suggests) to read the FASTA.
#'
#' @param ref_bed Path to a TSV BED table with a header, or a data frame, with
#'   columns `#chrom`, `start`, `end`, `target_name`, `length`, and `strand`.
#' @param target_fasta Path to a FASTA whose record names match `target_name`.
#' @param out Optional output TSV path. If `NULL`, results are returned without
#'   writing.
#' @param overwrite If `FALSE` (default), refuse to overwrite `out`.
#'
#' @return A tibble of `ref_bed` with a `ref_seq` column.
#'
#' @export
add_ref_seqs_with_targeted_ref_fasta <- function(ref_bed,
                                                 target_fasta,
                                                 out = NULL,
                                                 overwrite = FALSE) {
  options(dplyr.summarise.inform = FALSE)
  stop_if_output_exists(out, overwrite)
  check_suggested_pkg(
    "Biostrings",
    "reading targeted reference FASTA via add_ref_seqs_with_targeted_ref_fasta()"
  )

  if (is.character(ref_bed) && length(ref_bed) == 1L) {
    if (!file.exists(ref_bed)) {
      stop(ref_bed, " does not exist", call. = FALSE)
    }
    bed <- read_ref_bed_table(ref_bed)
    bed_label <- ref_bed
  } else if (is.data.frame(ref_bed)) {
    validate_ref_bed_table(ref_bed, "ref_bed")
    bed <- tibble::as_tibble(ref_bed)
    bed_label <- "ref_bed"
  } else {
    stop("`ref_bed` must be a file path or a data frame.", call. = FALSE)
  }

  stop_on_duplicate_names(bed$target_name, bed_label)

  dna <- Biostrings::readDNAStringSet(target_fasta)
  decomp <- set_decompose(bed$target_name, names(dna))
  if (length(decomp$only_in_vector_a) > 0) {
    stop(
      "the following loci were missing from the fasta file ", target_fasta,
      " but are in ", bed_label, "\n",
      paste(decomp$only_in_vector_a, collapse = ","),
      call. = FALSE
    )
  }

  dna_tab <- tibble::tibble(
    target_name = names(dna),
    ref_seq = unname(as.character(dna))
  )
  stop_on_duplicate_names(dna_tab$target_name, target_fasta)

  if ("ref_seq" %in% colnames(bed)) {
    bed <- dplyr::select(bed, -"ref_seq")
  }
  bed <- dplyr::left_join(bed, dna_tab, by = "target_name")

  if (!is.null(out)) {
    readr::write_tsv(bed, out)
  }
  bed
}
