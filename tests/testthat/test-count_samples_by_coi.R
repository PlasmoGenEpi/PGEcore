test_that("calculate_coi_distribution fills missing COI levels", {
  coi_calls <- tibble::tibble(
    specimen_name = c("s1", "s2", "s3"),
    coi = c(1, 3, 3)
  )
  out <- calculate_coi_distribution(coi_calls)
  expect_equal(out$coi, 1:3)
  expect_equal(out$n, c(1, 0, 2))
  expect_equal(out$proportion, c(1 / 3, 0, 2 / 3))
})

test_that("count_samples_by_coi works from data frame and file", {
  coi_calls <- tibble::tibble(
    specimen_name = c("s1", "s2"),
    coi = c(1.4, 2.6)
  )
  out <- count_samples_by_coi(coi_calls)
  expect_equal(out$coi, 1:3)
  expect_equal(out$n, c(1L, 0L, 1L))

  tmp_in <- tempfile(fileext = ".tsv")
  tmp_out <- tempfile(fileext = ".tsv")
  on.exit(unlink(c(tmp_in, tmp_out)), add = TRUE)
  readr::write_tsv(coi_calls, tmp_in)
  out2 <- count_samples_by_coi(tmp_in, output = tmp_out)
  expect_true(file.exists(tmp_out))
  expect_equal(out2$n, out$n)
})

test_that("count_samples_by_coi works on packaged example data", {
  path <- system.file("extdata", "example_coi_table.tsv", package = "PGEcore")
  skip_if_not(nzchar(path) && file.exists(path))
  out <- count_samples_by_coi(path)
  expect_true(all(c("coi", "n", "proportion") %in% names(out)))
  expect_equal(sum(out$n), 5)
})
