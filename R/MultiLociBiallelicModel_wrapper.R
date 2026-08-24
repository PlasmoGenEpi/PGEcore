#' Create Multi-Loci Biallelic Model input
#'
#' @param input_path Path to amino-acid calls TSV.
#' @param loci_group Path to loci-group TSV (`group_id`, `gene_id`, `aa_position`).
#' @param aa_sample_occurence_cut_off Amino-acid calls must occur in more than
#'   this number of samples to be included.
#' @return A list (`MLBM_object`) with `MLBM_data`, `staves_data`,
#'   `by_group_table`, and `groups`.
#' @keywords internal
create_MultiLociBiallelicModel_input <- function(input_path,
                                                 loci_group,
                                                 aa_sample_occurence_cut_off = 0) {
  print("Reading input data")
  original_input_data <- readr::read_tsv(
    input_path,
    col_types = readr::cols(specimen_name = readr::col_character())
  )

  if (aa_sample_occurence_cut_off > 0) {
    original_input_data <- original_input_data |>
      dplyr::group_by(.data$gene_id, .data$aa_position, .data$aa) |>
      dplyr::mutate(sample_count = dplyr::n_distinct(.data$specimen_name)) |>
      dplyr::filter(.data$sample_count > aa_sample_occurence_cut_off)
    input_data <- original_input_data |>
      dplyr::ungroup() |>
      dplyr::select(-"sample_count")
  } else {
    input_data <- original_input_data
  }

  MLBM_data <- input_data |>
    dplyr::mutate(identifier = paste(.data$gene_id, .data$aa_position, sep = ":")) |>
    dplyr::group_by(.data$specimen_name, .data$identifier) |>
    dplyr::mutate(value = dplyr::case_when(
      dplyr::n_distinct(.data$aa) > 2 ~ NA_real_,
      all(.data$ref_aa == .data$aa) ~ 0,
      all(.data$ref_aa != .data$aa) ~ 1,
      any(.data$ref_aa == .data$aa) & any(.data$ref_aa != .data$aa) ~ 2
    )) |>
    dplyr::ungroup() |>
    dplyr::select("specimen_name", "identifier", "value") |>
    dplyr::distinct() |>
    tidyr::pivot_wider(
      names_from = "identifier",
      values_from = "value",
      values_fill = NA
    ) |>
    dplyr::filter(!dplyr::if_any(dplyr::everything(), is.na))

  if (nrow(MLBM_data) == 0) {
    stop(
      "All samples have at least one missing genotype. Frequency estimation ",
      "is impossible.",
      call. = FALSE
    )
  }

  all_doubles <- all(vapply(MLBM_data[-1], is.double, logical(1)))
  if (!all_doubles) {
    stop(
      "Analysis object failed one or more validation checks:\n",
      "Amino acid table misformatted. Verify your inputs.",
      call. = FALSE
    )
  }

  staves_data <- input_data |>
    dplyr::mutate(
      staves = paste(.data$gene_id, .data$aa_position, .data$aa, sep = ":"),
      state = dplyr::if_else(
        substr(.data$staves, nchar(.data$staves), nchar(.data$staves)) == .data$ref_aa,
        0,
        1
      )
    ) |>
    dplyr::distinct(.data$staves, .keep_all = TRUE) |>
    dplyr::select("staves", "state") |>
    dplyr::mutate(prestave = sub(":([^:]+)$", "", .data$staves))

  MLBM_object <- list(
    MLBM_data = MLBM_data,
    staves_data = staves_data
  )
  loci_groups <- readr::read_tsv(
    loci_group,
    col_types = readr::cols(
      group_id = readr::col_character(),
      gene_id = readr::col_character(),
      aa_position = readr::col_integer()
    )
  )
  rules <- validate::validator(
    is.character(group_id),
    is.character(gene_id),
    is.integer(aa_position),
    !is.na(group_id),
    !is.na(gene_id),
    !is.na(aa_position)
  )
  print("Confronting input data with validation rules")
  stop_on_validate_fails(loci_groups, rules, "loci_group_table")

  identifiers_in_groups <- unique(
    paste0(loci_groups$gene_id, ":", loci_groups$aa_position)
  )
  input_data_aa_calls_check <- input_data |>
    dplyr::mutate(identifier = paste0(.data$gene_id, ":", .data$aa_position)) |>
    dplyr::filter(.data$identifier %in% identifiers_in_groups) |>
    dplyr::group_by(.data$gene_id, .data$aa_position) |>
    dplyr::summarise(
      distinct_aas = dplyr::n_distinct(.data$aa),
      aas = paste0(unique(sort(.data$aa)), collapse = ","),
      .groups = "drop"
    ) |>
    dplyr::filter(.data$distinct_aas > 2)

  if (nrow(input_data_aa_calls_check) > 0) {
    warning_message <- "The following positions are not biallelic "
    for (row in seq_len(nrow(input_data_aa_calls_check))) {
      warning_message <- paste0(
        warning_message,
        " ",
        input_data_aa_calls_check$gene_id[row],
        ":",
        input_data_aa_calls_check$aa_position[row],
        " aas: ",
        input_data_aa_calls_check$aas
      )
    }
    stop(warning_message, call. = FALSE)
  }

  match_group <- input_data |>
    dplyr::select("gene_id", "aa_position") |>
    dplyr::mutate(identifier = paste(.data$gene_id, .data$aa_position, sep = ":"))

  merged_data <- match_group |>
    dplyr::inner_join(
      loci_groups,
      by = c("gene_id", "aa_position"),
      relationship = "many-to-many"
    )

  grouped_list <- merged_data |>
    dplyr::group_by(.data$group_id) |>
    dplyr::summarise(
      identifiers = list(unique(.data$identifier)),
      .groups = "drop"
    )

  result_list <- stats::setNames(grouped_list$identifiers, grouped_list$group_id)
  result_tables <- lapply(names(result_list), function(group_name) {
    columns <- c("specimen_name", result_list[[group_name]])
    dplyr::select(MLBM_data, dplyr::all_of(columns))
  })
  names(result_tables) <- names(result_list)
  MLBM_object$by_group_table <- result_tables
  MLBM_object$groups <- unique(loci_groups$group_id)
  MLBM_object
}

#' Run the vendored Multi-Loci Biallelic Model MLE
#'
#' @param inputToMLE Matrix/data prepared for the vendored `mle()` function.
#' @return List with `lambda`, `plsf_table`, and `runtime`.
#' @keywords internal
run_MultiLociBiallelicModel <- function(inputToMLE) {
  runtime <- system.time({
    est <- .mlbm_vendor$mle(inputToMLE, id = TRUE)
  })
  plsf <- est$p
  plsf_table <- tibble::tibble(
    sequence = dimnames(plsf)[[2]],
    MLBM_frequency = c(plsf)
  )
  list(
    lambda = est$lambda,
    plsf_table = plsf_table,
    runtime = runtime
  )
}

#' Group variantstring-style mutations into stave format
#'
#' @param variant Character string of `;`-separated `gene:position:mutation`.
#' @return A stave-formatted string.
#' @keywords internal
make_stave <- function(variant) {
  parts <- stringr::str_split(variant, ";")[[1]]
  parsed <- parts |>
    purrr::map(~ stringr::str_match(.x, "(.*):(\\d+):(\\w+)")) |>
    purrr::map_dfr(~ tibble::tibble(gene = .x[2], position = .x[3], mutation = .x[4]))
  parsed |>
    dplyr::group_by(.data$gene) |>
    dplyr::summarize(
      positions = stringr::str_c(.data$position, collapse = "_"),
      mutations = stringr::str_c(.data$mutation, collapse = "_"),
      .groups = "drop"
    ) |>
    dplyr::mutate(
      formatted = stringr::str_c(.data$gene, ":", .data$positions, ":", .data$mutations)
    ) |>
    dplyr::pull(.data$formatted) |>
    stringr::str_c(collapse = ";")
}

#' Summarise MLBM haplotype frequencies as variant strings
#'
#' @param MLBM_res Result of `run_MultiLociBiallelicModel()`.
#' @param MLBM_object Result of `create_MultiLociBiallelicModel_input()`.
#' @param group_name Group id matching `MLBM_object$by_group_table`.
#' @return Tibble with `group_id`, `variant`, `freq`.
#' @keywords internal
summarise_MLBM_results <- function(MLBM_res, MLBM_object, group_name) {
  sequence_df <- MLBM_res$plsf_table |>
    dplyr::mutate(sequence_split = strsplit(.data$sequence, "")) |>
    tidyr::unnest_wider("sequence_split", names_sep = "_") |>
    dplyr::select(-"MLBM_frequency") |>
    dplyr::rename_with(
      ~ colnames(MLBM_object$by_group_table[[group_name]])[-1],
      dplyr::starts_with("sequence_split")
    )

  columns_to_replace <- colnames(sequence_df)[-1]
  sequence_long <- sequence_df |>
    dplyr::mutate(row_id = dplyr::row_number()) |>
    tidyr::pivot_longer(
      cols = dplyr::all_of(columns_to_replace),
      names_to = "prestave",
      values_to = "state"
    ) |>
    dplyr::mutate(state = as.integer(.data$state)) |>
    dplyr::left_join(
      MLBM_object$staves_data,
      by = c("prestave", "state"),
      relationship = "many-to-many"
    ) |>
    tidyr::separate("staves", into = c("gene", "pos", "aa"), sep = ":", convert = TRUE) |>
    dplyr::mutate(
      n_aa = 1,
      het = (.data$n_aa > 1),
      phased = FALSE,
      read_count = NA
    ) |>
    dplyr::arrange(.data$row_id, .data$gene, .data$pos)

  vs <- sequence_long |>
    dplyr::group_by(.data$row_id) |>
    tidyr::nest() |>
    dplyr::mutate(variant = purrr::map_chr(
      .data$data,
      ~ variantstring::long_to_variant(list(.x |>
        dplyr::select("gene", "pos", "n_aa", "het", "phased", "aa", "read_count") |>
        as.data.frame()))
    )) |>
    dplyr::select("row_id", "variant")

  MLBM_res$plsf_table <- MLBM_res$plsf_table |>
    dplyr::mutate(row_id = dplyr::row_number()) |>
    dplyr::left_join(vs, by = "row_id") |>
    dplyr::mutate(
      group_id = group_name,
      freq = .data$MLBM_frequency
    ) |>
    dplyr::select("group_id", "variant", "freq")

  MLBM_res$plsf_table
}

#' Estimate multilocus haplotype frequencies with MultiLociBiallelicModel
#'
#' File-oriented entry point used by the `MultiLociBiallelicModel_wrapper` CLI.
#' Requires suggested package **variantstring** 1.x. Samples with any missing
#' genotype are dropped (legacy behaviour; the model cannot handle missing data).
#'
#' @param aa_calls Path to amino-acid calls TSV (`specimen_name`, `gene_id`,
#'   `aa_position`, `ref_aa`, `aa`).
#' @param loci_group_table Path to loci groups TSV (`group_id`, `gene_id`,
#'   `aa_position`).
#' @param mlaf_output Output TSV (`group_id`, `variant`, `freq`).
#' @param aa_sample_occurence_cut_off Amino-acid calls must occur in more than
#'   this number of samples to be included (legacy default `0`).
#'
#' @return The bound MLAF tibble (invisibly after writing `mlaf_output`).
#' @export
MultiLociBiallelicModel_wrapper <- function(aa_calls,
                                            loci_group_table,
                                            mlaf_output,
                                            aa_sample_occurence_cut_off = 0) {
  if (is.null(aa_calls) || is.null(loci_group_table) || is.null(mlaf_output)) {
    stop(
      "Missing required arguments: --aa_calls, --loci_group_table, --mlaf_output",
      call. = FALSE
    )
  }
  check_variantstring_v1("STAVE strings via MultiLociBiallelicModel_wrapper()")

  MLBM_object <- create_MultiLociBiallelicModel_input(
    aa_calls,
    loci_group_table,
    aa_sample_occurence_cut_off
  )
  MLBM_res <- list()
  for (group_name in MLBM_object$groups) {
    print(paste0("Running MultiLociBiallelicModel on ", group_name))
    group_table <- MLBM_object$by_group_table[[group_name]]
    MLBM_tmp <- run_MultiLociBiallelicModel(group_table)
    MLBM_res[[group_name]] <- summarise_MLBM_results(
      MLBM_tmp,
      MLBM_object,
      group_name
    )
  }
  out <- dplyr::bind_rows(MLBM_res)
  readr::write_tsv(out, mlaf_output)
  invisible(out)
}
