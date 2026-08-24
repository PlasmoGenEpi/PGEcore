test_that("create_moire_input validates allele table columns", {
  tmp <- tempfile(fileext = ".tsv")
  on.exit(unlink(tmp), add = TRUE)
  writeLines("specimen_name\ttarget_name\nS1\tL1", tmp)

  err <- tryCatch(
    create_moire_input(
      tmp, TRUE, 10, 10, FALSE, 1, 1, 1, 1, 1, 1, 0.1, 10, 2, 2,
      FALSE, 1, 0, 1, TRUE, Inf
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
    tmp, FALSE, 2, 2, FALSE, 1, 1, 1, 1, 1, 1, 0.1, 10, 2, 2,
    FALSE, 1, 0, 1, TRUE, 1
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
    dcifer_ibd_wrapper(allele_table = NULL, btwn_host_rel_output = "x.tsv"),
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
    btwn_host_rel_output = out,
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
      loci_groups_input = "x",
      mlaf_output = "y",
      coi_output = "z"
    ),
    "missing the following arguments|snp.slicer|variantstring"
  )
})

test_that("snpslice_wrapper integration skips without snp.slicer", {
  skip_if_not_installed("snp.slicer")
  skip_if_not_installed("variantstring")

  aa <- system.file("extdata", "example2_amino_acid_calls.tsv", package = "PGEcore")
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
      groups = "g",
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

  aa <- system.file("extdata", "example_amino_acid_calls.tsv", package = "PGEcore")
  groups <- system.file("extdata", "example_loci_groups.tsv", package = "PGEcore")
  skip_if_not(nzchar(aa) && file.exists(aa))
  skip_if_not(nzchar(groups) && file.exists(groups))
  tbl <- read_fem_aa_calls(aa)
  expect_true(all(c("specimen_name", "gene_id", "aa") %in% names(tbl)))
})
