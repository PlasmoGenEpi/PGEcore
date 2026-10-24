#' Load amino acid calls for naive allele prevalence estimation
#'
#' @param path Path to a TSV of amino acid calls.
#' @return A tibble with columns `specimen_name`, `target_name`, and `variant`.
#' @keywords internal
prev_parse_aa_calls <- function(path) {
  readr::read_tsv(
    path,
    col_types = readr::cols(
      specimen_name = readr::col_character(),
      gene_id = readr::col_character(),
      reads = readr::col_integer(),
      aa_position = readr::col_integer(),
      ref_aa = readr::col_character(),
      aa = readr::col_character()
    ),
    col_select = c(
      "specimen_name", "gene_id",
      "aa_position", "aa"
    )
  ) |>
    tidyr::unite("target_name", "gene_id", "aa_position", sep = ":") |>
    dplyr::rename(variant = "aa")
}

#' Load microhaplotype calls for naive allele prevalence estimation
#'
#' @param path Path to a TSV of microhaplotype genotypes.
#' @return A tibble with columns including `specimen_name`, `target_name`, and
#'   `variant`.
#' @keywords internal
prev_parse_mh_calls <- function(path) {
  readr::read_tsv(
    path,
    col_types = readr::cols(
      specimen_name = readr::col_character(),
      target_name = readr::col_character(),
      seq = readr::col_character(),
      reads = readr::col_integer()
    )
  ) |>
    dplyr::rename(variant = "seq")
}

#' Calculate allele prevalence
#'
#' Presence-based allele prevalence: the fraction of specimens at a locus in
#' which each allele is detected. A mixed infection contributes to every allele
#' it carries, so prevalences at a locus can sum to more than 1.
#'
#' @param allele_table Allele table with `specimen_name`, `target_name`, and
#'   `variant`.
#' @return A tibble with `target_name`, `variant`, `prev`, `sample_count`, and
#'   `sample_total`.
#' @keywords internal
calculate_prevalence <- function(allele_table) {
  allele_table |>
    dplyr::group_by(.data$target_name) |>
    dplyr::mutate(sample_total = dplyr::n_distinct(.data$specimen_name)) |>
    dplyr::group_by(.data$target_name, .data$variant, .data$sample_total) |>
    dplyr::summarise(
      sample_count = dplyr::n_distinct(.data$specimen_name)
    ) |>
    dplyr::mutate(prev = .data$sample_count / .data$sample_total) |>
    dplyr::relocate("prev", .before = "sample_total")
}

#' Collapse alleles in a specimen into one observed genotype
#'
#' Alleles are sorted and joined with `/` (malaria surveillance notation for
#' mixed calls, e.g. `K/T`).
#'
#' @param allele_table Allele table with `specimen_name`, `target_name`, and
#'   `variant`.
#' @return A tibble with one row per specimen and locus.
#' @keywords internal
collapse_to_observed_genotype <- function(allele_table) {
  allele_table |>
    dplyr::distinct(.data$specimen_name, .data$target_name, .data$variant) |>
    dplyr::group_by(.data$specimen_name, .data$target_name) |>
    dplyr::summarise(
      variant = paste(sort(.data$variant), collapse = "/"),
      .groups = "drop"
    )
}

#' Calculate observed-genotype prevalence
#'
#' Each specimen contributes to exactly one genotype: the set of alleles
#' detected at that locus. Prevalences at a locus sum to 1.
#'
#' @param allele_table Allele table with `specimen_name`, `target_name`, and
#'   `variant`.
#' @return A tibble with `target_name`, `variant`, `prev`, `sample_count`, and
#'   `sample_total`.
#' @keywords internal
calculate_genotype_prevalence <- function(allele_table) {
  calculate_prevalence(collapse_to_observed_genotype(allele_table))
}

#' Format naive allele-prevalence output for AA or MH input
#'
#' @param prevalence Prevalence table with `target_name` and `variant`.
#' @param from_aa Logical; `TRUE` when input was amino acid calls.
#' @return Formatted tibble for writing.
#' @keywords internal
format_naive_prev_output <- function(prevalence, from_aa) {
  if (isTRUE(from_aa)) {
    prevalence |>
      tidyr::separate_wider_delim(
        "target_name",
        ":",
        names = c("gene_id", "aa_position")
      ) |>
      dplyr::rename(aa = "variant") |>
      convert_single_locus_table_to_stave(
        additional_columns = c("prev", "sample_count", "sample_total")
      )
  } else {
    dplyr::rename(prevalence, seq = "variant")
  }
}

#' Estimate allele or observed-genotype prevalence from AA or microhaplotype calls
#'
#' Exactly one of `aa_calls` or `allele_table` must be provided. The two
#' methods count mixed infections differently.
#'
#' **`allele`** (default) is *allele prevalence*: the fraction of specimens in
#' which each allele is detected. A mixed infection counts toward **every**
#' allele it carries, so prevalences at a locus can sum to more than 1
#' (Hastings and Smith 2008; *Malaria Journal* 7:130).
#'
#' **`genotype`** is *observed genotype prevalence*: each specimen counts
#' **once**, as the full set of alleles detected at that locus. Mixed
#' infections are written with `/` (`A/B`, `K/T`, `CVMNK/CVIET`). Prevalences
#' at a locus sum to 1.
#'
#' Example — three specimens at one locus (monoclonal A, mixed A and B,
#' monoclonal B):
#'
#' | `method` | A | A/B | B | sum |
#' | -------- | - | --- | - | --- |
#' | `allele` | 2/3 | — | 2/3 | 4/3 |
#' | `genotype` | 1/3 | 1/3 | 1/3 | 1 |
#'
#' ## Inputs
#'
#' - **`aa_calls`** (optional): AA calls (`specimen_name`, `gene_id`,
#'   `aa_position`, `aa`). See `vignette("input-formats", package = "PGEcore")`.
#' - **`allele_table`** (optional): Allele table (`specimen_name`,
#'   `target_name`, `seq`). See the same vignette.
#'
#' ## Outputs
#'
#' - **`output`** (optional): Prevalence TSV. For AA input: `variant`, `prev`,
#'   `sample_count`, `sample_total` (`variant` is STAVE
#'   `gene_id:aa_position:aa`; mixed genotypes use `aa` like `C/R`). For
#'   microhaplotype input: `target_name`, `seq`, `prev`, `sample_count`,
#'   `sample_total` (mixed genotypes join sequences with `/`). If `NULL`,
#'   results are returned without writing a file.
#'
#' ## Running
#'
#' ```r
#' estimate_allele_prevalence_naive(
#'   aa_calls = "aa_calls.tsv",
#'   output = "prevalence.tsv"
#' )
#' estimate_allele_prevalence_naive(
#'   aa_calls = "aa_calls.tsv",
#'   method = "genotype",
#'   output = "genotype_prevalence.tsv"
#' )
#' ```
#'
#' ```bash
#' Rscript exec/estimate_allele_prevalence_naive \
#'   --aa_calls aa_calls.tsv \
#'   --output prevalence.tsv
#' Rscript exec/estimate_allele_prevalence_naive \
#'   --aa_calls aa_calls.tsv \
#'   --method genotype \
#'   --output genotype_prevalence.tsv
#' ```
#'
#' @param aa_calls Optional path to an AA calls TSV. See *Inputs*.
#' @param allele_table Optional path to an allele table TSV. See *Inputs*.
#' @param method `"allele"` (default): each mixed infection counts toward every
#'   allele it carries (allele prevalence). `"genotype"`: each specimen counts
#'   once as its observed allele combination (observed genotype prevalence).
#' @param output Optional output TSV path. Default for the CLI is
#'   `prevalence.tsv`.
#'
#' @return A tibble of prevalences. For amino acid input: `variant`,
#'   `prev`, `sample_count`, `sample_total`. For microhaplotype input:
#'   `target_name`, `seq`, `prev`, `sample_count`, `sample_total`.
#'
#' @seealso `vignette("input-formats", package = "PGEcore")`
#'
#' @references
#' Hastings IM, Smith TA (2008). MalHaploFreq: A computer programme for
#' estimating malaria haplotype frequencies from blood samples.
#' *Malaria Journal* 7: 130.
#'
#' @examples
#' aa_path <- system.file(
#'   "extdata", "example_aa_calls.tsv",
#'   package = "PGEcore"
#' )
#' estimate_allele_prevalence_naive(aa_calls = aa_path)
#' estimate_allele_prevalence_naive(aa_calls = aa_path, method = "genotype")
#'
#' @export
estimate_allele_prevalence_naive <- function(aa_calls = NULL,
                                             allele_table = NULL,
                                             method = "allele",
                                             output = NULL) {
  options(dplyr.summarise.inform = FALSE)

  if (is.null(aa_calls) == is.null(allele_table)) {
    stop(
      "One and only one of the args --aa_calls and --allele_table should be ",
      "provided.",
      call. = FALSE
    )
  }
  if (!method %in% c("allele", "genotype")) {
    stop(method, " is not a valid method", call. = FALSE)
  }

  from_aa <- !is.null(aa_calls)
  if (from_aa) {
    if (!file.exists(aa_calls)) {
      stop(aa_calls, " does not exist", call. = FALSE)
    }
    allele_table <- prev_parse_aa_calls(aa_calls)
  } else {
    if (!file.exists(allele_table)) {
      stop(allele_table, " does not exist", call. = FALSE)
    }
    allele_table <- prev_parse_mh_calls(allele_table)
  }

  prevalence <- switch(
    method,
    allele = calculate_prevalence(allele_table),
    genotype = calculate_genotype_prevalence(allele_table)
  )
  prev_output <- format_naive_prev_output(prevalence, from_aa = from_aa)

  if (!is.null(output)) {
    readr::write_tsv(prev_output, output)
  }

  prev_output
}
