test_that("multilocus_prevfreq_naive_variantstring checks suggested package", {
  err <- tryCatch(
    multilocus_prevfreq_naive_variantstring(
      aa_table = "missing.tsv",
      loci_groups_input = "also_missing.tsv"
    ),
    error = function(e) conditionMessage(e)
  )
  expect_true(grepl("variantstring|does not exist", err))
})

test_that("multilocus_prevfreq_naive_variantstring works on packaged example data", {
  skip_if_not_installed("variantstring")

  aa_path <- system.file(
    "extdata",
    "example_amino_acid_calls.tsv",
    package = "PGEcore"
  )
  groups_path <- system.file(
    "extdata",
    "example_loci_groups.tsv",
    package = "PGEcore"
  )
  skip_if_not(nzchar(aa_path) && file.exists(aa_path))
  skip_if_not(nzchar(groups_path) && file.exists(groups_path))

  out <- multilocus_prevfreq_naive_variantstring(aa_path, groups_path)
  expect_equal(
    names(out),
    c("group_id", "variant", "prev", "freq", "sample_total")
  )
  expect_true(nrow(out) > 0)
  expect_true(all(out$prev >= 0 & out$prev <= 1))
  expect_true(all(out$freq >= 0 & out$freq <= 1))

  tmp_out <- tempfile(fileext = ".tsv")
  on.exit(unlink(tmp_out), add = TRUE)
  written <- multilocus_prevfreq_naive_variantstring(
    aa_path,
    groups_path,
    output_path = tmp_out
  )
  expect_true(file.exists(tmp_out))
  expect_equal(nrow(written), nrow(out))
})
