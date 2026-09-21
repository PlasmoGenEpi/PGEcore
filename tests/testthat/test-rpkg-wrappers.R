test_that("create_moire_input validates allele table columns", {
  tmp <- tempfile(fileext = ".tsv")
  on.exit(unlink(tmp), add = TRUE)
  writeLines("specimen_name\ttarget_name\nS1\tL1", tmp)

  err <- tryCatch(
    create_moire_input(
      tmp, TRUE, 10, 10, 1, FALSE, 1, 1, 1, 1, 1, 1, 0.1, 10, 2, 2,
      FALSE, 1, 0, 1, TRUE, Inf, 1, 1
    ),
    error = function(e) conditionMessage(e)
  )
  expect_true(grepl("checkmate|moire|Missing required columns|seq", err))
})

test_that("run_moire / moire_wrapper skip when moire is unavailable", {
  skip_if_not_installed("moire")
  skip_if_not_installed("checkmate")

  path <- system.file("extdata", "example2_allele_table.tsv", package = "PGEcore")
  skip_if_not(nzchar(path) && file.exists(path))

  dat <- readr::read_tsv(path, show_col_types = FALSE)
  dat <- dat[dat$specimen_name %in% head(unique(dat$specimen_name), 2), ]
  dat <- dat[dat$target_name %in% head(unique(dat$target_name), 3), ]
  tmp <- tempfile(fileext = ".tsv")
  on.exit(unlink(tmp), add = TRUE)
  readr::write_tsv(dat, tmp)

  obj <- create_moire_input(
    tmp, FALSE, 2, 2, 1, FALSE, 1, 1, 1, 1, 1, 1, 0.1, 10, 2, 2,
    FALSE, 1, 0, 1, TRUE, 1, 1, 1
  )
  expect_true(all(c("moire_data", "moire_parameters") %in% names(obj)))
})

test_that("malariaem_wrapper validates missing allele table", {
  expect_error(
    malariaem_wrapper(allele_table = tempfile()),
    "does not exist|malaria.em|checkmate"
  )
})

test_that("run_malariaem integration skips when malaria.em is unavailable", {
  skip_if_not_installed("malaria.em")
  skip_if_not_installed("checkmate")

  path <- system.file("extdata", "example_target_groups.tsv", package = "PGEcore")
  skip_if_not(nzchar(path) && file.exists(path))
  groups <- load_malariaem_target_groups(path)
  expect_true(all(c("group_id", "target_name") %in% names(groups)))
})

test_that("dcifer_slaf_wrapper validates required paths", {
  expect_error(
    dcifer_slaf_wrapper(allele_table = NULL, slaf_output = "x.tsv"),
    "required|dcifer"
  )
})

test_that("dcifer_slaf_wrapper integration skips without dcifer", {
  skip_if_not_installed("dcifer")
  path <- system.file("extdata", "example_allele_table.tsv", package = "PGEcore")
  skip_if_not(nzchar(path) && file.exists(path))
  out <- tempfile(fileext = ".tsv")
  on.exit(unlink(out), add = TRUE)
  slaf <- dcifer_slaf_wrapper(allele_table = path, slaf_output = out)
  expect_true(file.exists(out))
  expect_true("freq" %in% names(slaf))
})

test_that("dcifer_ibd_wrapper validates required paths", {
  expect_error(
    dcifer_ibd_wrapper(allele_table = NULL, relatedness_output = "x.tsv"),
    "required|dcifer|foreach"
  )
})

test_that("dcifer_ibd_wrapper integration skips without dcifer", {
  skip_if_not_installed("dcifer")
  skip_if_not_installed("foreach")
  skip_if_not_installed("doParallel")
  skip_if_not_installed("parallelly")
  skip_if_not_installed("iterators")

  path <- system.file("extdata", "example_allele_table.tsv", package = "PGEcore")
  skip_if_not(nzchar(path) && file.exists(path))
  dat <- readr::read_tsv(path, show_col_types = FALSE)
  dat <- dat[
    dat$specimen_name %in% head(unique(dat$specimen_name), 3) &
      dat$target_name %in% head(unique(dat$target_name), 4),
  ]
  tmp <- tempfile(fileext = ".tsv")
  out <- tempfile(fileext = ".tsv")
  on.exit(unlink(c(tmp, out)), add = TRUE)
  readr::write_tsv(dat, tmp)
  res <- dcifer_ibd_wrapper(
    allele_table = tmp,
    relatedness_output = out,
    threads = 1L
  )
  expect_true(file.exists(out))
  expect_true("btwn_host_rel" %in% names(
    dplyr::rename(res, btwn_host_rel = "estimate")
  ) || "estimate" %in% names(res))
})

test_that("snpslice_wrapper validates required arguments", {
  expect_error(
    snpslice_wrapper(
      allele_table = NULL,
      loci_groups = "x",
      mlaf_output = "y",
      coi_output = "z"
    ),
    "missing the following arguments|snp.slicer|variantstring"
  )
})

test_that("snpslice_wrapper rejects an unknown estimator", {
  expect_error(
    snpslice_wrapper(
      allele_table = "a",
      loci_groups = "x",
      mlaf_output = "y",
      coi_output = "z",
      estimator = "mcmc"
    ),
    "estimator must be one of|snp.slicer|variantstring"
  )
})

test_that("snpslice_wrapper integration skips without snp.slicer", {
  skip_if_not_installed("snp.slicer")
  skip_if_not_installed("variantstring")

  aa <- system.file("extdata", "example2_aa_calls.tsv", package = "PGEcore")
  lg <- system.file("extdata", "example_loci_groups.tsv", package = "PGEcore")
  skip_if_not(nzchar(aa) && file.exists(aa))
  skip_if_not(nzchar(lg) && file.exists(lg))

  expect_error(
    create_snpslice_allele_table_input(
      aa,
      target_name_col = "aa_locus",
      target_value_col = "aa",
      target_count_col = "reads"
    ),
    NA
  )
})

test_that("FreqEstimationModel_wrapper validates required arguments", {
  expect_error(
    FreqEstimationModel_wrapper(
      aa_calls = NULL,
      coi = 1,
      loci_groups = "g",
      mlaf_output = "o"
    ),
    "required|FreqEstimationModel|variantstring"
  )
})

test_that("calculate_avg_COI accepts numeric COI or a table", {
  expect_equal(calculate_avg_COI("2"), 2)
  path <- system.file("extdata", "example_coi_table.tsv", package = "PGEcore")
  skip_if_not(nzchar(path) && file.exists(path))
  avg <- calculate_avg_COI(path)
  expect_true(is.numeric(avg) && length(avg) == 1)
})

test_that("FreqEstimationModel_wrapper integration skips without FEM", {
  skip_if_not_installed("FreqEstimationModel")
  skip_if_not_installed("variantstring")
  skip_if_not_installed("plyr")
  skip_if_not_installed("coda")
  skip_if_not_installed("abind")
  skip_if_not_installed("foreach")
  skip_if_not_installed("doMC")

  aa <- system.file("extdata", "example_aa_calls.tsv", package = "PGEcore")
  groups <- system.file("extdata", "example_loci_groups.tsv", package = "PGEcore")
  skip_if_not(nzchar(aa) && file.exists(aa))
  skip_if_not(nzchar(groups) && file.exists(groups))
  tbl <- read_fem_aa_calls(aa)
  expect_true(all(c("specimen_name", "gene_id", "aa") %in% names(tbl)))
})

# The multinomial observation model keeps targets with more than two alleles;
# every other model drops them from the loci groups. The example data has a
# triallelic dhps 540 (E/K/N) inside the pfdhps and pfdhfr_pfdhps groups.
test_that("snpslice_wrapper keeps multi-allelic targets only for the multinomial model", {
  skip_if_not_installed("snp.slicer")
  skip_if_not_installed("variantstring")
  skip_if_not(
    exists("snp_slice_multinomial", asNamespace("snp.slicer")),
    "installed snp.slicer lacks the multinomial model"
  )
  aa <- system.file("extdata", "example2_aa_calls.tsv", package = "PGEcore")
  lg <- system.file("extdata", "example_loci_groups.tsv", package = "PGEcore")
  skip_if_not(nzchar(aa) && file.exists(aa))
  skip_if_not(nzchar(lg) && file.exists(lg))
  out <- withr::local_tempdir()

  allele_tbl <- create_snpslice_allele_table_input(
    aa, target_name_col = "aa_locus", target_value_col = "aa", target_count_col = "reads"
  )
  # One warning per group that contains the triallelic codon.
  w <- testthat::capture_warnings(
    groups_bi <- create_snpslice_loci_group_input(lg, allele_tbl, target_name_col = "aa_locus")
  )
  expect_true(all(grepl("more than two alleles", w)))
  expect_length(w, 2L)
  expect_false("PF3D7_0810800.1:540" %in% groups_bi$pfdhps)
  expect_false("PF3D7_0810800.1:540" %in% groups_bi$pfdhfr_pfdhps)
  expect_no_warning(
    groups_multi <- create_snpslice_loci_group_input(
      lg, allele_tbl, target_name_col = "aa_locus", allow_multiallelic = TRUE
    )
  )
  expect_true("PF3D7_0810800.1:540" %in% groups_multi$pfdhps)

  res <- snpslice_wrapper(
    aa, lg,
    mlaf_output = file.path(out, "mlaf.tsv"),
    coi_output = file.path(out, "coi.tsv"),
    convergence_output = file.path(out, "conv.tsv"),
    model = "multinomial", loci_limit = 12, n_sample = 30, n_burnin = 30,
    n_chains = 1, seed = 3, estimator = "map"
  )
  dhps <- res$mlaf[res$mlaf$group_id == "pfdhps", ]
  alleles_540 <- sub(".*540:[A-Z]_([A-Z])$", "\\1", dhps$variant)
  # A third allele at 540 survives into the multilocus haplotypes.
  expect_true(all(c("E", "K", "N") %in% alleles_540))
  expect_equal(sum(dhps$freq), 1, tolerance = 1e-8)
  expect_true(file.exists(file.path(out, "mlaf.tsv")))
})

test_that("snpslice_wrapper rejects an unknown model", {
  expect_error(
    snpslice_wrapper(
      allele_table = "a", loci_groups = "x", mlaf_output = "y", coi_output = "z",
      model = "dirichlet"
    ),
    "model must be one of|snp.slicer|variantstring"
  )
})

# Ju et al. (2024) count identical haplotypes once. snpslice_dedup_matrices()
# merges strains with identical dictionary rows; a specimen carrying any of the
# duplicates carries the merged strain once.
test_that("snpslice_dedup_matrices merges identical dictionary rows", {
  A <- rbind(c(1, 1, 0), c(0, 1, 1), c(1, 0, 0))
  D <- rbind(c(0, 1), c(0, 1), c(1, 0))   # strains 1 and 2 identical
  d <- snpslice_dedup_matrices(A, D)
  expect_equal(d$D, rbind(c(0, 1), c(1, 0)))
  expect_equal(unname(d$A), rbind(c(1, 0), c(1, 1), c(1, 0)))
  expect_equal(rowSums(d$A), c(1, 2, 1))       # specimen 1: two copies -> one
  # Nothing to merge: unchanged.
  D2 <- rbind(c(0, 1), c(1, 1), c(1, 0))
  d2 <- snpslice_dedup_matrices(A, D2)
  expect_equal(d2$D, D2)
  expect_equal(unname(d2$A), A)
  # Empty dictionary passes through.
  d0 <- snpslice_dedup_matrices(A[, 0, drop = FALSE], D[0, , drop = FALSE])
  expect_equal(ncol(d0$A), 0)
})

test_that("snpslice_wrapper reports coi_chain_mean and averaged_freq, and dedup_haplotypes is honoured", {
  skip_if_not_installed("snp.slicer")
  skip_if_not_installed("variantstring")
  skip_if_not(
    exists("snp_slice_multinomial", asNamespace("snp.slicer")),
    "installed snp.slicer lacks the multinomial model"
  )
  aa <- system.file("extdata", "example2_aa_calls.tsv", package = "PGEcore")
  lg <- system.file("extdata", "example_loci_groups.tsv", package = "PGEcore")
  skip_if_not(nzchar(aa) && file.exists(aa))
  out <- withr::local_tempdir()
  run <- function(dedup) {
    snpslice_wrapper(
      aa, lg,
      mlaf_output = file.path(out, "mlaf.tsv"), coi_output = file.path(out, "coi.tsv"),
      convergence_output = file.path(out, "conv.tsv"),
      model = "multinomial", loci_limit = 12, n_sample = 30, n_burnin = 30, gap = 30,
      n_chains = 2, seed = 3, estimator = "final_sample", dedup_haplotypes = dedup
    )
  }
  res <- run(FALSE)
  expect_equal(names(res$coi), c("specimen_name", "coi", "coi_cons_weighted", "coi_chain_mean"))
  expect_true(all(c("freq", "averaged_freq", "allele_count", "allele_total") %in% names(res$mlaf)))
  sums <- res$mlaf |> dplyr::group_by(.data$group_id) |>
    dplyr::summarise(f = sum(.data$freq), a = sum(.data$averaged_freq), .groups = "drop")
  expect_equal(sums$f, rep(1, nrow(sums)), tolerance = 1e-8)
  expect_equal(sums$a, rep(1, nrow(sums)), tolerance = 1e-8)
  expect_true(all(res$mlaf$freq > 0 | res$mlaf$averaged_freq > 0))
  # The chain mean is bounded by the per-chain values, and with 2 chains is a
  # multiple of 0.5 for the final-sample estimator.
  expect_true(all(res$coi$coi_chain_mean >= 1))
  expect_true(all(abs(res$coi$coi_chain_mean * 2 - round(res$coi$coi_chain_mean * 2)) < 1e-8))

  res_dd <- run(TRUE)
  expect_equal(names(res_dd$coi), names(res$coi))
  # De-duplication can only lower strain counts, never raise them.
  expect_true(all(res_dd$coi$coi <= res$coi$coi))
  expect_true(all(res_dd$mlaf$allele_total <= max(res$mlaf$allele_total)))
  # Dedup is the default.
  expect_true(formals(snpslice_wrapper)$dedup_haplotypes)
})

