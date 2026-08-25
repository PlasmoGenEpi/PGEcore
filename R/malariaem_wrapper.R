#' Load and format an allele table for malaria.em
#'
#' @keywords internal
load_malariaem_allele_table <- function(file_path) {
  if (!file.exists(file_path)) {
    stop(file_path, " does not exist", call. = FALSE)
  }

  message("Reading input data...")
  mhap <- readr::read_tsv(
    file_path,
    col_types = readr::cols(
      specimen_name = readr::col_character(),
      target_name = readr::col_character(),
      seq = readr::col_character()
    ),
    col_select = c("specimen_name", "target_name", "seq")
  )

  rules <- validate::validator(
    is.character(specimen_name),
    is.character(target_name),
    is.character(seq),
    !is.na(specimen_name),
    !is.na(target_name),
    !is.na(seq)
  )
  stop_on_validate_fails(mhap, rules, "allele_table")

  message("Converting to matrix...")
  mhap |>
    dplyr::select("specimen_name", "target_name", "seq") |>
    dplyr::group_by(.data$specimen_name, .data$target_name) |>
    tidyr::pivot_wider(
      names_from = "target_name",
      values_from = "seq",
      values_fn = ~ paste(unique(.), collapse = " ")
    ) |>
    tibble::column_to_rownames(var = "specimen_name") |>
    as.matrix()
}

#' Load target groups for subsetted malaria.em runs
#'
#' @keywords internal
load_malariaem_target_groups <- function(file_path) {
  if (!file.exists(file_path)) {
    stop(file_path, " does not exist", call. = FALSE)
  }

  message("Reading target groups...")
  target_groups <- readr::read_tsv(
    file_path,
    col_types = readr::cols(
      group_id = readr::col_character(),
      target_name = readr::col_character()
    ),
    show_col_types = FALSE
  )

  rules <- validate::validator(
    is.character(group_id),
    is.character(target_name),
    !is.na(group_id),
    !is.na(target_name)
  )
  stop_on_validate_fails(target_groups, rules, "target_groups")
  target_groups
}

#' Run malaria.em on one allele matrix
#'
#' @keywords internal
run_malariaem_helper <- function(matrix,
                                 test_size,
                                 max_size,
                                 label = NULL) {
  min_targets <- max(apply(matrix, 2, malaria.em::num.unique.allele.locus))
  if (identical(test_size, "min")) {
    coi_range <- seq_len(min_targets)
  } else {
    top <- as.numeric(test_size)
    if (min_targets > top) {
      top <- min_targets
    }
    coi_range <- seq_len(top)
  }

  max_size_n <- as.numeric(max_size)
  if (max(coi_range) > max_size_n) {
    stop("Too many alleles to run EM efficiently", call. = FALSE)
  }

  matrix <- matrix[rowSums(is.na(matrix)) == 0, , drop = FALSE]

  sample_name <- matrix |>
    as.data.frame() |>
    tibble::rownames_to_column("specimen_name") |>
    dplyr::distinct(.data$specimen_name) |>
    tibble::rowid_to_column("ids")

  target_names <- colnames(matrix)

  locus_allele_counts <- apply(matrix, 2, function(col) {
    alleles <- unlist(strsplit(col, "\\s+"))
    alleles <- alleles[!is.na(alleles) & nzchar(alleles)]
    length(unique(alleles))
  })

  if (all(locus_allele_counts == 1)) {
    message(
      "All loci", if (!is.null(label)) paste0(" [", label, "]"),
      " have only a single microhaplotype across the population; there is ",
      "only one possible haplotype. Skipping malaria.em and phasing every ",
      "sample as the only possible haplotype."
    )

    single_alleles <- vapply(target_names, function(tn) {
      alleles <- unlist(strsplit(matrix[, tn], "\\s+"))
      alleles <- alleles[!is.na(alleles) & nzchar(alleles)]
      unique(alleles)[1]
    }, character(1))

    single_hap <- tibble::tibble(
      target_name = target_names,
      seq = unname(single_alleles)
    )

    gt_freq_summary <- single_hap |>
      dplyr::mutate(gt_id = 1L, freq = 1, freq_se = 0) |>
      dplyr::select("gt_id", "target_name", "seq", "freq", "freq_se")

    gt_phase_summary <- sample_name |>
      dplyr::select("specimen_name") |>
      tidyr::crossing(single_hap) |>
      dplyr::mutate(gt_id = 1L, posterior_est = 1, phase_id = 1L) |>
      dplyr::select(
        "specimen_name", "target_name", "seq", "gt_id", "posterior_est", "phase_id"
      )

    if (!is.null(label)) {
      gt_freq_summary <- gt_freq_summary |> dplyr::mutate(group_id = label)
      gt_phase_summary <- gt_phase_summary |> dplyr::mutate(group_id = label)
    }

    return(list(
      output = NULL,
      gt_freq_summary = gt_freq_summary,
      gt_phase_summary = gt_phase_summary,
      group_id = label
    ))
  }

  message("Running malaria.em", if (!is.null(label)) paste0(" [", label, "]"), "...")
  output <- malaria.em::malaria.em(
    matrix,
    sizes = coi_range,
    locus.label = target_names
  )

  message("Summarizing population-level multi-locus genotype frequency estimates + SE...")
  gt_freq_summary <- output$haplo.prob.tab |>
    as.data.frame() |>
    tibble::rowid_to_column("gt_id") |>
    tidyr::pivot_longer(
      cols = -c("gt_id", "hap.prob", "hap.prob.std"),
      names_to = "hap_id",
      values_to = "seq"
    ) |>
    dplyr::select(
      "gt_id",
      target_name = "hap_id",
      "seq",
      freq = "hap.prob",
      freq_se = "hap.prob.std"
    ) |>
    dplyr::mutate(
      freq = as.numeric(.data$freq),
      freq_se = as.numeric(.data$freq_se)
    )

  if (!is.null(label)) {
    gt_freq_summary <- gt_freq_summary |> dplyr::mutate(group_id = label)
  }

  message("Summarizing sample-level phased multi-locus genotypes + posterior probability estimates...")
  haplo_pred <- as.data.frame(output$pred.haplo.set)
  if ("output$pred.haplo.set" %in% names(haplo_pred)) {
    haplo_pred <- dplyr::rename(haplo_pred, haplo.set = `output$pred.haplo.set`)
  } else if (!"haplo.set" %in% names(haplo_pred) && ncol(haplo_pred) == 1L) {
    names(haplo_pred) <- "haplo.set"
  }
  haplo_pred <- tibble::rowid_to_column(haplo_pred, "ids")
  haplo_pred_prob <- as.data.frame(output$haplo.sets)

  gt_phase_summary <- haplo_pred |>
    as.data.frame() |>
    dplyr::left_join(haplo_pred_prob, by = c("ids", "haplo.set")) |>
    dplyr::left_join(sample_name, by = "ids") |>
    dplyr::select("specimen_name", "haplo.set", posterior_est = "post.p") |>
    tidyr::separate_rows("haplo.set", sep = " ") |>
    dplyr::rename(gt_id = "haplo.set") |>
    dplyr::mutate(gt_id = as.integer(.data$gt_id)) |>
    dplyr::left_join(
      gt_freq_summary,
      by = "gt_id",
      relationship = "many-to-many"
    ) |>
    dplyr::select(
      "specimen_name", "target_name", "seq", "gt_id", "posterior_est"
    ) |>
    dplyr::mutate(phase_id = dplyr::dense_rank(.data$gt_id), .by = "specimen_name")

  if (!is.null(label)) {
    gt_phase_summary <- gt_phase_summary |> dplyr::mutate(group_id = label)
  }

  list(
    output = output,
    gt_freq_summary = gt_freq_summary,
    gt_phase_summary = gt_phase_summary,
    group_id = label
  )
}

#' Run malaria.em and write frequency and phase summaries
#'
#' The **malaria.em** package is an optional dependency (Suggests). It is not
#' installed automatically with PGEcore.
#'
#' @param matrix Allele matrix (specimens as rows, loci as columns).
#' @param test_size Maximum COI to test, or `"min"` for the minimum allowed.
#' @param max_size Error if inferred COI range exceeds this cutoff.
#' @param subset_targets If `TRUE`, run separately for each `group_id`.
#' @param target_groups Data frame with `group_id` and `target_name`.
#' @param freq_output Path for genotype-frequency TSV.
#' @param phase_output Path for phasing TSV.
#'
#' @return A list of malaria.em results (or a named list per group).
#' @export
run_malariaem <- function(matrix,
                          test_size = "min",
                          max_size = "8",
                          subset_targets = FALSE,
                          target_groups = NULL,
                          freq_output = NULL,
                          phase_output = NULL) {
  check_suggested_pkg("malaria.em", "EM haplotype inference via run_malariaem()")
  check_suggested_pkg("checkmate", "malaria.em input validation")

  checkmate::assert_matrix(matrix)
  checkmate::assert_flag(subset_targets)
  if (subset_targets) {
    checkmate::assert_data_frame(target_groups, min.rows = 1, min.cols = 2)
    checkmate::assert_subset(c("group_id", "target_name"), choices = names(target_groups))
  }

  freq_outdir <- if (is.null(freq_output)) getwd() else dirname(freq_output)
  phase_outdir <- if (is.null(phase_output)) getwd() else dirname(phase_output)
  if (!dir.exists(freq_outdir)) {
    dir.create(freq_outdir, recursive = TRUE, showWarnings = FALSE)
  }
  if (!dir.exists(phase_outdir)) {
    dir.create(phase_outdir, recursive = TRUE, showWarnings = FALSE)
  }
  checkmate::assert_directory_exists(freq_outdir, access = "w")
  checkmate::assert_directory_exists(phase_outdir, access = "w")

  if (!subset_targets) {
    res <- run_malariaem_helper(matrix, test_size, max_size, label = NULL)
    readr::write_tsv(res$gt_freq_summary, freq_output)
    readr::write_tsv(res$gt_phase_summary, phase_output)
    return(res)
  }

  group_ids <- unique(target_groups$group_id) |> as.character()
  results_list <- stats::setNames(vector("list", length(group_ids)), group_ids)

  for (gid in group_ids) {
    targets <- unique(target_groups$target_name[target_groups$group_id == gid])
    if (length(targets) == 0L) {
      warning("No targets found for group_id '", gid, "'. Skipping.", call. = FALSE)
      next
    }
    non_existing <- setdiff(targets, colnames(matrix))
    if (length(non_existing) > 0) {
      stop(
        "Group '", gid, "' requests targets not present in matrix: ",
        paste(non_existing, collapse = ", "),
        call. = FALSE
      )
    }
    matrix_subset <- matrix[, targets, drop = FALSE]
    if (ncol(matrix_subset) == 0) {
      warning("Matrix for group_id '", gid, "' is empty; skipping.", call. = FALSE)
      next
    }
    results_list[[gid]] <- run_malariaem_helper(
      matrix_subset,
      test_size,
      max_size,
      label = gid
    )
  }

  all_gt_freq_summary <- dplyr::bind_rows(
    lapply(results_list, `[[`, "gt_freq_summary")
  )
  all_gt_phase_summary <- dplyr::bind_rows(
    lapply(results_list, `[[`, "gt_phase_summary")
  )
  readr::write_tsv(all_gt_freq_summary, freq_output)
  readr::write_tsv(all_gt_phase_summary, phase_output)
  results_list
}

#' Run malaria.em from allele-table and output paths
#'
#' File-oriented entry point used by the `malariaem_wrapper` CLI.
#'
#' @param allele_table Path to allele TSV (`specimen_name`, `target_name`, `seq`).
#' @param subset_targets Logical; subset by `target_groups`.
#' @param target_groups Optional path to groups TSV.
#' @param max_size COI cutoff (legacy CLI default `"8"`).
#' @param test_size COI size to test, or `"min"`.
#' @param freq_output Frequency summary path.
#' @param phase_output Phase summary path.
#' @param seed Random seed.
#'
#' @return The object returned by [run_malariaem()].
#' @export
malariaem_wrapper <- function(allele_table,
                              subset_targets = FALSE,
                              target_groups = NULL,
                              max_size = "8",
                              test_size = "min",
                              freq_output = "gt_freq_summary_all.tsv",
                              phase_output = "gt_phase_summary_all.tsv",
                              seed = 1L) {
  options(dplyr.summarise.inform = FALSE)
  check_suggested_pkg("malaria.em", "malaria.em via malariaem_wrapper()")
  check_suggested_pkg("checkmate", "malaria.em input validation")

  set.seed(as.numeric(seed))
  matrix <- load_malariaem_allele_table(allele_table)

  tg <- NULL
  if (isTRUE(subset_targets)) {
    if (is.null(target_groups)) {
      stop("--target_groups is required when --subset_targets is TRUE", call. = FALSE)
    }
    tg <- load_malariaem_target_groups(target_groups)
  }

  run_malariaem(
    matrix = matrix,
    test_size = test_size,
    max_size = max_size,
    subset_targets = isTRUE(subset_targets),
    target_groups = tg,
    freq_output = freq_output,
    phase_output = phase_output
  )
}
