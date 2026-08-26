#' Read amino acid calls for variantstring multilocus prev/freq
#'
#' @param aa_calls Path to a TSV of amino acid calls.
#' @return A tibble with `specimen_name`, `gene`, `pos`, `reads`, `aa`, `n_aa`.
#' @keywords internal
read_mlp_vs_aa_table <- function(aa_calls) {
  stopifnot(is.character(aa_calls), length(aa_calls) == 1L)
  if (!file.exists(aa_calls)) {
    stop(aa_calls, " does not exist", call. = FALSE)
  }

  df_aa <- utils::read.table(
    aa_calls,
    header = TRUE,
    colClasses = c(specimen_name = "character")
  )

  validate_required_columns(
    df_aa,
    c("specimen_name", "gene_id", "aa_position", "reads", "aa"),
    "aa_calls"
  )

  rules <- validate::validator(
    is.character(specimen_name),
    is.character(gene_id),
    is.integer(aa_position),
    is.integer(reads),
    is.character(aa),
    !is.na(specimen_name),
    !is.na(gene_id),
    !is.na(aa_position),
    !is.na(reads),
    !is.na(aa)
  )
  stop_on_validate_fails(df_aa, rules, "aa_calls")

  df_aa |>
    dplyr::select("specimen_name", "gene_id", "aa_position", "reads", "aa") |>
    dplyr::rename(gene = "gene_id", pos = "aa_position") |>
    dplyr::group_by(.data$specimen_name, .data$gene, .data$pos) |>
    dplyr::mutate(n_aa = dplyr::n()) |>
    dplyr::ungroup()
}

#' Read loci-group definitions for variantstring multilocus prev/freq
#'
#' @param loci_groups_path Path to a TSV with `group_id`, `gene_id`, `aa_position`.
#' @return A tibble of loci groups.
#' @keywords internal
read_mlp_vs_loci_groups <- function(loci_groups_path) {
  stopifnot(is.character(loci_groups_path), length(loci_groups_path) == 1L)
  if (!file.exists(loci_groups_path)) {
    stop(loci_groups_path, " does not exist", call. = FALSE)
  }

  loci_groups <- readr::read_tsv(
    loci_groups_path,
    col_types = readr::cols(
      .default = readr::col_character(),
      aa_position = readr::col_integer()
    ),
    progress = FALSE
  )

  validate_required_columns(
    loci_groups,
    c("group_id", "gene_id", "aa_position"),
    "loci_groups"
  )

  rules <- validate::validator(
    is.character(group_id),
    is.character(gene_id),
    is.integer(aa_position),
    !is.na(group_id),
    !is.na(gene_id),
    !is.na(aa_position)
  )
  stop_on_validate_fails(loci_groups, rules, "loci_groups")
  loci_groups
}

#' Convert long-form amino acid calls to variant strings
#'
#' @param aa_table Tibble with `specimen_name`, `gene`, `pos`, `reads`, `aa`, `n_aa`.
#' @return Character vector of variant strings, one per specimen.
#' @keywords internal
aa_table_to_variant <- function(aa_table) {
  long_list <- aa_table |>
    dplyr::mutate(
      het = (.data$n_aa > 1),
      phased = FALSE
    ) |>
    dplyr::select(
      "gene", "pos", "n_aa", "het", "phased", "aa",
      read_count = "reads"
    ) |>
    split(f = aa_table$specimen_name)
  names(long_list) <- NULL

  variantstring::long_to_variant(long_list)
}

#' Extract unambiguous component genotypes from variant strings
#'
#' @param variant_strings Character vector of variant strings.
#' @return Unique component variant strings with no heterozygous loci.
#' @keywords internal
extract_component_variants <- function(variant_strings) {
  variantstring::check_variant_string(variant_strings)

  pos_unique <- unique(variantstring::position_from_variant_string(variant_strings))

  vs_allpos <- mapply(
    function(x) {
      stats::na.omit(
        unique(
          variantstring::subset_position(
            position_string = x,
            variant_strings = variant_strings
          )
        )
      )
    },
    pos_unique,
    SIMPLIFY = FALSE
  ) |>
    unlist() |>
    unique()

  stats::na.omit(unique(unlist(variantstring::get_component_variants(vs_allpos))))
}

#' Prevalence and frequency of target variants vs a dataset
#'
#' @param target_variants Target variant strings.
#' @param comparison_variants Dataset variant strings.
#' @return A data frame with `variant`, `prev`, `freq`, `sample_total`.
#' @keywords internal
calculate_variant_prevalence <- function(target_variants, comparison_variants) {
  variantstring::check_variant_string(target_variants)
  variantstring::check_variant_string(comparison_variants)

  l <- list()
  for (i in seq_along(target_variants)) {
    df_match <- variantstring::compare_variant_string(
      target_string = target_variants[[i]],
      comparison_strings = comparison_variants
    )
    df_match$match_pos <- variantstring::compare_position_string(
      target_string = variantstring::position_from_variant_string(
        target_variants[[i]]
      ),
      comparison_strings = comparison_variants
    )

    l[[i]] <- df_match |>
      dplyr::filter(.data$match_pos) |>
      dplyr::filter(!.data$ambiguous) |>
      dplyr::summarise(
        variant = target_variants[[i]],
        prev = mean(.data$match),
        freq = mean(.data$prop),
        sample_total = dplyr::n()
      )
  }

  dplyr::bind_rows(l)
}

#' Compute prev/freq for one loci group via variantstring
#'
#' @param group_loci Tibble with `gene_id` and `aa_position`.
#' @param aa_table Amino acid calls in long form.
#' @return Tibble with `variant`, `prev`, `freq`, `sample_total`.
#' @keywords internal
compute_prevfreq_for_group <- function(group_loci, aa_table) {
  variant_strings <- aa_table |>
    dplyr::inner_join(
      group_loci,
      by = c("gene" = "gene_id", "pos" = "aa_position")
    ) |>
    aa_table_to_variant()

  calculate_variant_prevalence(
    extract_component_variants(variant_strings),
    variant_strings
  )
}

#' Write a variantstring prev/freq table
#'
#' @param df_prev Data frame with `group_id`, `variant`, `prev`, `freq`, `sample_total`.
#' @param output Output TSV path.
#' @keywords internal
write_mlp_vs_prev <- function(df_prev, output) {
  stopifnot(is.data.frame(df_prev))
  stopifnot(
    all(names(df_prev) == c("group_id", "variant", "prev", "freq", "sample_total"))
  )
  stopifnot(is.character(output), length(output) == 1L)
  readr::write_tsv(df_prev, file = output)
}

#' Estimate multilocus prevalence and frequency with variantstring
#'
#' Converts amino acid calls to
#' [variant strings](https://github.com/mrc-ide/variantstring), extracts
#' unambiguous component genotypes, and estimates prevalence and frequency of
#' each component against the full dataset.
#'
#' ## Inputs
#'
#' - **`aa_calls`**: AA calls (`specimen_name`, `gene_id`, `aa_position`,
#'   `reads`, `aa`). See `vignette("input-formats", package = "PGEcore")`.
#' - **`loci_groups`**: Loci groups (`group_id`, `gene_id`, `aa_position`).
#'   See the same vignette.
#'
#' ## Outputs
#'
#' - **`output`** (optional): Prev/freq TSV with columns `group_id`,
#'   `variant`, `prev`, `freq`, and `sample_total`. If `NULL`, results are
#'   returned without writing.
#'
#' ## Running
#'
#' ```r
#' multilocus_prevfreq_naive_variantstring(
#'   aa_calls = "aa_calls.tsv",
#'   loci_groups = "loci_groups.tsv",
#'   output = "multilocus_prevfreq.tsv"
#' )
#' ```
#'
#' ```bash
#' Rscript exec/multilocus_prevfreq_naive_variantstring \
#'   --aa_calls aa_calls.tsv \
#'   --loci_groups loci_groups.tsv \
#'   --output multilocus_prevfreq.tsv
#' ```
#'
#' Requires the optional **variantstring** package (Suggests). It is not
#' installed automatically with PGEcore.
#'
#' @param aa_calls Path to an AA calls TSV. See *Inputs*.
#' @param loci_groups Path to a loci groups TSV. See *Inputs*.
#' @param output Optional path for the prev/freq TSV.
#'
#' @return A tibble with columns `group_id`, `variant`, `prev`, `freq`, and
#'   `sample_total`.
#'
#' @seealso `vignette("input-formats", package = "PGEcore")`
#'
#' @examplesIf requireNamespace("variantstring", quietly = TRUE)
#' aa_path <- system.file(
#'   "extdata", "example_aa_calls.tsv",
#'   package = "PGEcore"
#' )
#' groups_path <- system.file(
#'   "extdata", "example_loci_groups.tsv",
#'   package = "PGEcore"
#' )
#' multilocus_prevfreq_naive_variantstring(aa_path, groups_path)
#'
#' @export
multilocus_prevfreq_naive_variantstring <- function(aa_calls,
                                                    loci_groups,
                                                    output = NULL) {
  options(dplyr.summarise.inform = FALSE)
  check_suggested_pkg(
    "variantstring",
    "multilocus prev/freq via multilocus_prevfreq_naive_variantstring()"
  )

  aa_calls <- read_mlp_vs_aa_table(aa_calls)
  loci_groups <- read_mlp_vs_loci_groups(loci_groups)

  df_prev <- loci_groups |>
    dplyr::select("group_id", "gene_id", "aa_position") |>
    tidyr::nest(group_loci = c("gene_id", "aa_position")) |>
    dplyr::mutate(
      group_mlafp = lapply(
        .data$group_loci,
        compute_prevfreq_for_group,
        aa_calls
      )
    ) |>
    dplyr::select(-"group_loci") |>
    tidyr::unnest("group_mlafp")

  if (!is.null(output)) {
    write_mlp_vs_prev(df_prev = df_prev, output = output)
  }

  df_prev
}
