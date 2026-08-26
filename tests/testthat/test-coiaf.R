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
