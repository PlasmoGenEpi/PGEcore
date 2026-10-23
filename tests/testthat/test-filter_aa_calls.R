example_aa_calls <- function() {
  tibble::tibble(
    specimen_name = "s1",
    gene_id = "G1",
    aa_position = c(1, 1, 1, 2),
    aa = c("N", "R", "I", "C"),
    reads = c(80, 15, 5, 50)
  )
}

test_that("filter_aa_calls keeps calls that pass min_reads", {
  out <- filter_aa_calls(example_aa_calls(), min_reads = 20)
  expect_equal(out$aa, c("N", "C"))
  expect_equal(out$reads, c(80, 50))
})

test_that("filter_aa_calls keeps calls that pass min_wsaf", {
  out <- filter_aa_calls(example_aa_calls(), min_wsaf = 0.10)
  expect_equal(out$aa, c("N", "R", "C"))
})

test_that("filter_aa_calls requires calls to pass both thresholds", {
  out <- filter_aa_calls(
    example_aa_calls(),
    min_reads = 20,
    min_wsaf = 0.10
  )
  expect_equal(out$aa, c("N", "C"))
})

test_that("filter_aa_calls applies thresholds independently (no WSAF renormalize)", {
  # R WSAF is 15/100 = 0.15. After dropping I (reads < 6), a sequential
  # renormalize would give R 15/95 ≈ 0.158 and keep it at min_wsaf = 0.151.
  out <- filter_aa_calls(
    example_aa_calls(),
    min_reads = 6,
    min_wsaf = 0.151
  )
  expect_equal(out$aa, c("N", "C"))
})

test_that("filter_aa_calls is a no-op when no thresholds are set", {
  calls <- example_aa_calls()
  out <- filter_aa_calls(calls)
  expect_equal(out, calls)
})

test_that("filter_aa_calls keeps extra columns", {
  calls <- dplyr::mutate(example_aa_calls(), gene = "g", target_name = "t1")
  out <- filter_aa_calls(calls, min_reads = 20)
  expect_true(all(c("gene", "target_name") %in% names(out)))
  expect_false("wsaf" %in% names(out))
})

test_that("filter_aa_calls validates columns, thresholds, and overwrite", {
  expect_error(
    filter_aa_calls(data.frame(specimen_name = "s")),
    "Missing required columns"
  )
  expect_error(
    filter_aa_calls(example_aa_calls(), min_reads = -1),
    "min_reads"
  )
  expect_error(
    filter_aa_calls(example_aa_calls(), min_wsaf = 1.2),
    "min_wsaf"
  )
  expect_error(
    filter_aa_calls(example_aa_calls(), min_reads = 1000),
    "All amino-acid calls were removed"
  )

  tmp <- tempfile(fileext = ".tsv")
  on.exit(unlink(tmp), add = TRUE)
  readr::write_tsv(example_aa_calls(), tmp)
  expect_error(
    filter_aa_calls(example_aa_calls(), output = tmp, overwrite = FALSE),
    "already exists"
  )
})

test_that("filter_aa_calls works on packaged example", {
  path <- system.file(
    "extdata",
    "example_aa_calls.tsv",
    package = "PGEcore"
  )
  skip_if_not(nzchar(path) && file.exists(path))
  out <- filter_aa_calls(path, min_reads = 3)
  expect_true(nrow(out) >= 1)
  expect_true(all(out$reads >= 3))
  expect_true(all(
    c("specimen_name", "gene_id", "aa_position", "aa", "reads") %in% names(out)
  ))
})
