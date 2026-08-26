test_that("calculate_af_presence_absence and read_count_prop agree on toy data", {
  allele_table <- tibble::tibble(
    specimen_name = c("s1", "s1", "s2"),
    target_name = c("L1", "L1", "L1"),
    variant = c("A", "B", "A"),
    reads = c(10L, 10L, 5L)
  )

  pa <- calculate_af_presence_absence(allele_table)
  expect_equal(pa$freq[pa$variant == "A"], 2 / 3)
  expect_equal(pa$freq[pa$variant == "B"], 1 / 3)

  rc <- calculate_af_read_count_prop(allele_table)
  expect_equal(rc$freq[rc$variant == "A"], 0.75)
  expect_equal(rc$freq[rc$variant == "B"], 0.25)
})

test_that("estimate_allele_frequency_naive validates exclusive inputs and method", {
  expect_error(
    estimate_allele_frequency_naive(),
    "One and only one"
  )
  expect_error(
    estimate_allele_frequency_naive(
      aa_calls = "missing.tsv",
      allele_table = "also_missing.tsv"
    ),
    "One and only one"
  )

  aa_path <- system.file(
    "extdata",
    "example_aa_calls.tsv",
    package = "PGEcore"
  )
  skip_if_not(nzchar(aa_path) && file.exists(aa_path))
  expect_error(
    estimate_allele_frequency_naive(aa_calls = aa_path, method = "nope"),
    "not a valid method"
  )
})

test_that("estimate_allele_frequency_naive works on packaged example data", {
  aa_path <- system.file(
    "extdata",
    "example_aa_calls.tsv",
    package = "PGEcore"
  )
  mh_path <- system.file(
    "extdata",
    "example_allele_table.tsv",
    package = "PGEcore"
  )
  skip_if_not(nzchar(aa_path) && file.exists(aa_path))
  skip_if_not(nzchar(mh_path) && file.exists(mh_path))

  aa_out <- estimate_allele_frequency_naive(
    aa_calls = aa_path,
    method = "presence_absence"
  )
  expect_true(all(c("variant", "freq") %in% names(aa_out)))
  expect_true(all(grepl(":", aa_out$variant)))

  mh_out <- estimate_allele_frequency_naive(
    allele_table = mh_path,
    method = "read_count_prop"
  )
  expect_true(all(c("target_name", "seq", "freq") %in% names(mh_out)))

  tmp_out <- tempfile(fileext = ".tsv")
  on.exit(unlink(tmp_out), add = TRUE)
  written <- estimate_allele_frequency_naive(
    aa_calls = aa_path,
    output = tmp_out
  )
  expect_true(file.exists(tmp_out))
  expect_equal(nrow(written), nrow(aa_out))
})
