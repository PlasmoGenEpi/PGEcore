test_that("run_coiaf validates inputs without requiring coiaf for arg checks", {
  # Parameter validation happens after check_suggested_pkg.
  # If coiaf is missing, we get the Suggests error; if present, invalid args error.
  bad <- data.frame(specimen_name = "s1")
  err <- tryCatch(
    run_coiaf(bad, seq_error = 2),
    error = function(e) conditionMessage(e)
  )
  expect_true(
    grepl("coiaf|seq_error|Missing required columns", err)
  )
})

test_that("run_coiaf integration skips when coiaf is unavailable", {
  skip_if_not_installed("coiaf")

  snp_path <- system.file(
    "extdata",
    "example_collapsed_snp_calls.tsv",
    package = "PGEcore"
  )
  skip_if_not(nzchar(snp_path) && file.exists(snp_path))

  snp_calls <- readr::read_tsv(snp_path, show_col_types = FALSE, n_max = 2000)
  specimens <- head(unique(snp_calls$specimen_name), 2)
  snp_calls <- snp_calls[snp_calls$specimen_name %in% specimens, ]

  expect_error(run_coiaf(snp_calls, seq_error = -1), "seq_error")

  # Auto-calculate PLMAF from the same subset to avoid mismatch warnings
  result <- suppressWarnings(run_coiaf(snp_calls, plmaf = NULL, max_coi = 5))
  expect_true(all(c("specimen_name", "coi_freq", "coi_variant") %in% names(result)))
  expect_true(nrow(result) >= 1)
})

test_that("coiaf_optimize resolves the monoclonal sentinel from its attribute", {
  # coiaf signals "no variant loci, COI is 1" by returning NaN with the answer in
  # an estimated_coi attribute rather than as the value.
  local_mocked_bindings(
    optimize_coi = function(...) {
      structure(NaN, notes = "Too few variant loci", estimated_coi = 1,
                num_variant_loci = 0)
    },
    .package = "coiaf"
  )
  expect_equal(coiaf_optimize(data.frame()), 1)
})

test_that("coiaf_optimize keeps the fitted value when the sentinel is not 1", {
  # The second sentinel path flags a low variant count but still carries the
  # optimiser's own estimate, which must not be collapsed to 1.
  local_mocked_bindings(
    optimize_coi = function(...) {
      structure(NaN, notes = "Too few variant loci", estimated_coi = 2.75,
                num_variant_loci = 3, expected_num_loci = 9)
    },
    .package = "coiaf"
  )
  expect_equal(coiaf_optimize(data.frame()), 2.75)
})

test_that("coiaf_optimize leaves an unexplained NaN as NaN", {
  # A NaN with no attribute is a real failure and must not become a number.
  local_mocked_bindings(optimize_coi = function(...) NaN, .package = "coiaf")
  expect_true(is.nan(coiaf_optimize(data.frame())))
})

test_that("coiaf_optimize passes ordinary estimates through", {
  local_mocked_bindings(optimize_coi = function(...) 3.5, .package = "coiaf")
  expect_equal(coiaf_optimize(data.frame()), 3.5)
})

test_that("run_coiaf weights loci by depth, not by the minor-allele count", {
  skip_if_not_installed("coiaf")
  # coiaf documents `coverage` as the read depth at each locus and uses it as the
  # per-locus weight, so a homozygous locus must still carry its full depth.
  seen <- NULL
  local_mocked_bindings(
    optimize_coi = function(data, ...) {
      seen <<- rbind(seen, as.data.frame(data))
      1
    },
    .package = "coiaf"
  )
  snp_calls <- data.frame(
    specimen_name = c("s1", "s1", "s1", "s1"),
    snp_name      = c("L1", "L1", "L2", "L2"),
    seq_base      = c("A", "T", "C", "G"),
    reads         = c(90L, 10L, 100L, 0L)
  )
  invisible(run_coiaf(snp_calls))
  # L1 depth 100 (minor allele T, 10 reads); L2 depth 100 (minor allele G, 0 reads)
  expect_true(all(seen$coverage == 100))
  expect_false(any(seen$coverage %in% c(10, 0)))
})
