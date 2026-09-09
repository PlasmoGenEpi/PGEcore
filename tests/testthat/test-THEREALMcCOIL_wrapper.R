test_that("THEREALMcCOIL_wrapper requires input and output paths", {
  expect_error(
    THEREALMcCOIL_wrapper(slaf_output = "s.tsv", coi_output = "c.tsv"),
    "snp_calls"
  )
  expect_error(
    THEREALMcCOIL_wrapper(
      snp_calls = "in.tsv",
      slaf_output = NULL,
      coi_output = "c.tsv"
    ),
    "slaf_output"
  )
  expect_error(
    THEREALMcCOIL_wrapper(
      snp_calls = "in.tsv",
      slaf_output = "s.tsv",
      coi_output = NULL
    ),
    "coi_output"
  )
})

test_that("read_and_preprocess_snp_call validates columns", {
  tmp <- tempfile(fileext = ".tsv")
  on.exit(unlink(tmp), add = TRUE)
  writeLines("specimen_name\tsnp_name\nS1\tL1", tmp)
  expect_error(
    PGEcore:::read_and_preprocess_snp_call(tmp),
    "Missing required columns|reads|seq_base"
  )
})

test_that("prep_input_categorical recodes major/minor/het/missing", {
  df <- tibble::tibble(
    specimen_name = c("S1", "S1", "S2", "S2", "S3"),
    snp_name = c("L1", "L1", "L1", "L2", "L1"),
    seq_base = c("A", "T", "A", "G", "T"),
    reads = c(10L, 8L, 12L, 5L, 9L)
  )
  mat <- PGEcore:::prep_input_categorical(df)
  expect_true(is.data.frame(mat))
  expect_setequal(rownames(mat), c("S1", "S2", "S3"))
  expect_true("L1" %in% colnames(mat))
  expect_true(all(as.matrix(mat) %in% c(-1, 0, 0.5, 1)))
  expect_equal(unname(mat["S1", "L1"]), 0.5)
  expect_equal(unname(mat["S2", "L1"]), 1)
  expect_equal(unname(mat["S3", "L1"]), 0)
  expect_equal(unname(mat["S2", "L2"]), 1)
  expect_equal(unname(mat["S1", "L2"]), -1)
})

test_that("prep_input_prop returns paired allele count matrices", {
  df <- tibble::tibble(
    specimen_name = c("S1", "S1", "S2", "S2"),
    snp_name = c("L1", "L1", "L1", "L1"),
    seq_base = c("A", "T", "A", "T"),
    reads = c(10L, 2L, 3L, 7L)
  )
  mats <- PGEcore:::prep_input_prop(df)
  expect_named(mats, c("a1", "a2"))
  expect_equal(nrow(mats$a1), 2)
  expect_equal(ncol(mats$a1), 1)
  expect_equal(mats$a1 + mats$a2, mats$a1 + mats$a2)
  expect_equal(as.numeric(mats$a1["S1", "L1"] + mats$a2["S1", "L1"]), 12)
})

test_that("prep_input_prop drops non-biallelic loci and indexes within locus", {
  df <- tibble::tibble(
    specimen_name = c("S1", "S1", "S1", "S1", "S2", "S2", "S2"),
    snp_name = c("L1", "L2", "L2", "L3", "L1", "L2", "L3"),
    seq_base = c("A", "A", "T", "C", "A", "T", "G"),
    reads = c(5L, 10L, 2L, 4L, 6L, 8L, 3L)
  )
  expect_warning(
    mats <- PGEcore:::prep_input_prop(df),
    "1 monomorphic"
  )
  expect_setequal(colnames(mats$a1), c("L2", "L3"))
  # Allele 1 is the first allele of that locus alphabetically, independent of
  # how many alleles the other loci have.
  expect_equal(as.numeric(mats$a1["S1", "L2"]), 10)
  expect_equal(as.numeric(mats$a2["S1", "L2"]), 2)
  expect_equal(as.numeric(mats$a1["S1", "L3"]), 4)
  expect_equal(as.numeric(mats$a2["S2", "L3"]), 3)
})

test_that("filter_biallelic errors when no biallelic loci remain", {
  df <- tibble::tibble(
    specimen_name = c("S1", "S2"),
    snp_name = c("L1", "L1"),
    seq_base = c("A", "A"),
    reads = c(5L, 6L)
  )
  expect_error(
    suppressWarnings(PGEcore:::filter_biallelic(df)),
    "No biallelic loci"
  )
})

test_that("run_mccoil_categorical stops when n or k is too small", {
  small <- as.data.frame(matrix(c(0, 1, 0.5, 1), nrow = 2, ncol = 2))
  rownames(small) <- c("S1", "S2")
  colnames(small) <- c("L1", "L2")
  expect_error(
    PGEcore:::run_mccoil_categorical(
      small,
      threshold_ind = 1,
      threshold_site = 1,
      totalrun = 5,
      burnin = 1,
      path = tempdir()
    ),
    "Sample size is too small"
  )
})

test_that("THEREALMcCOIL_wrapper short MCMC on example data", {
  skip_on_cran()
  skip_if_not_installed("posterior")
  path <- system.file(
    "extdata",
    "example_collapsed_snp_calls.tsv",
    package = "PGEcore"
  )
  skip_if_not(nzchar(path) && file.exists(path))
  df <- PGEcore:::read_and_preprocess_snp_call(path)
  cat_in <- PGEcore:::prep_input_categorical(df)
  skip_if(
    nrow(cat_in) <= 10 || ncol(cat_in) <= 10,
    "example SNP table is too small for McCOIL (needs n>10 and k>10)"
  )
  slaf <- tempfile(fileext = ".tsv")
  coi <- tempfile(fileext = ".tsv")
  convergence <- tempfile(fileext = ".tsv")
  on.exit(unlink(c(slaf, coi, convergence)), add = TRUE)
  # The toy run keeps too few draws for stable ESS estimates, so posterior
  # warns about capping them.
  res <- suppressWarnings(THEREALMcCOIL_wrapper(
    snp_calls = path,
    slaf_output = slaf,
    coi_output = coi,
    convergence_output = convergence,
    model = "categorical",
    maxCOI = 5L,
    threshold_ind = 5L,
    threshold_site = 5L,
    totalrun = 30L,
    burnin = 5L,
    M0 = 2L,
    n_chains = 2L
  ))
  expect_true(file.exists(slaf))
  expect_true(file.exists(coi))
  expect_true(file.exists(convergence))
  slaf_df <- readr::read_tsv(slaf, show_col_types = FALSE)
  coi_df <- readr::read_tsv(coi, show_col_types = FALSE)
  convergence_df <- readr::read_tsv(convergence, show_col_types = FALSE)
  expect_true(all(c("variant", "freq") %in% names(slaf_df)))
  expect_true(all(c("specimen_name", "coi") %in% names(coi_df)))
  expect_true(all(c("variable", "rhat", "ess_bulk") %in% names(convergence_df)))
  expect_gt(nrow(slaf_df), 0)
  expect_gt(nrow(coi_df), 0)
  expect_equal(nrow(convergence_df), nrow(slaf_df) + nrow(coi_df))
  expect_true("slaf" %in% names(res))
})

test_that("THEREALMcCOIL_wrapper is reproducible for a given seed", {
  skip_on_cran()
  skip_if_not_installed("posterior")
  path <- system.file(
    "extdata",
    "example_collapsed_snp_calls.tsv",
    package = "PGEcore"
  )
  skip_if_not(nzchar(path) && file.exists(path))
  df <- PGEcore:::read_and_preprocess_snp_call(path)
  cat_in <- PGEcore:::prep_input_categorical(df)
  skip_if(
    nrow(cat_in) <= 10 || ncol(cat_in) <= 10,
    "example SNP table is too small for McCOIL (needs n>10 and k>10)"
  )
  run_once <- function(seed) {
    slaf <- tempfile(fileext = ".tsv")
    coi <- tempfile(fileext = ".tsv")
    convergence <- tempfile(fileext = ".tsv")
    on.exit(unlink(c(slaf, coi, convergence)), add = TRUE)
    suppressWarnings(THEREALMcCOIL_wrapper(
      snp_calls = path,
      slaf_output = slaf,
      coi_output = coi,
      convergence_output = convergence,
      model = "categorical",
      maxCOI = 5L,
      threshold_ind = 5L,
      threshold_site = 5L,
      totalrun = 30L,
      burnin = 5L,
      M0 = 2L,
      n_chains = 1L,
      seed = seed
    ))$coi
  }
  expect_equal(run_once(321L), run_once(321L))
})
