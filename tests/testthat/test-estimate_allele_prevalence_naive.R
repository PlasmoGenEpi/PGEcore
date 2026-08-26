test_that("calculate_prevalence counts specimens per allele", {
  allele_table <- tibble::tibble(
    specimen_name = c("s1", "s1", "s2"),
    target_name = c("L1", "L1", "L1"),
    variant = c("A", "B", "A")
  )
  prev <- calculate_prevalence(allele_table)
  expect_equal(prev$prev[prev$variant == "A"], 1)
  expect_equal(prev$prev[prev$variant == "B"], 0.5)
  expect_equal(prev$sample_total[1], 2)
})

test_that("estimate_allele_prevalence_naive validates exclusive inputs", {
  expect_error(
    estimate_allele_prevalence_naive(),
    "One and only one"
  )
  expect_error(
    estimate_allele_prevalence_naive(
      aa_calls = "missing.tsv",
      allele_table = "also_missing.tsv"
    ),
    "One and only one"
  )
})

test_that("estimate_allele_prevalence_naive works on packaged example data", {
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

  aa_out <- estimate_allele_prevalence_naive(aa_calls = aa_path)
  expect_true(
    all(c("variant", "prev", "sample_count", "sample_total") %in% names(aa_out))
  )
  expect_false("gene_id" %in% names(aa_out))

  mh_out <- estimate_allele_prevalence_naive(allele_table = mh_path)
  expect_true(
    all(
      c("target_name", "seq", "prev", "sample_count", "sample_total") %in%
        names(mh_out)
    )
  )

  tmp_out <- tempfile(fileext = ".tsv")
  on.exit(unlink(tmp_out), add = TRUE)
  written <- estimate_allele_prevalence_naive(
    aa_calls = aa_path,
    output = tmp_out
  )
  expect_true(file.exists(tmp_out))
  expect_equal(nrow(written), nrow(aa_out))
})
