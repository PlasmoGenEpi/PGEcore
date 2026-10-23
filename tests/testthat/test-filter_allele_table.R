example_calls <- function() {
  tibble::tibble(
    specimen_name = "s1",
    target_name = c("t1", "t1", "t1", "t2"),
    seq = c("AAA", "TTT", "CCC", "GGG"),
    reads = c(80, 15, 5, 50)
  )
}

test_that("filter_allele_table keeps alleles that pass min_reads", {
  out <- filter_allele_table(example_calls(), min_reads = 20)
  expect_equal(out$seq, c("AAA", "GGG"))
  expect_equal(out$reads, c(80, 50))
})

test_that("filter_allele_table keeps alleles that pass min_wsaf", {
  out <- filter_allele_table(example_calls(), min_wsaf = 0.10)
  expect_equal(out$seq, c("AAA", "TTT", "GGG"))
})

test_that("filter_allele_table requires alleles to pass both thresholds", {
  out <- filter_allele_table(
    example_calls(),
    min_reads = 20,
    min_wsaf = 0.10
  )
  expect_equal(out$seq, c("AAA", "GGG"))
})

test_that("filter_allele_table applies thresholds independently (no WSAF renormalize)", {
  # TTT WSAF is 15/100 = 0.15. After dropping CCC (reads < 6), a sequential
  # renormalize would give TTT 15/95 ≈ 0.158 and keep it at min_wsaf = 0.151.
  out <- filter_allele_table(
    example_calls(),
    min_reads = 6,
    min_wsaf = 0.151
  )
  expect_equal(out$seq, c("AAA", "GGG"))
})

test_that("filter_allele_table is a no-op when no thresholds are set", {
  calls <- example_calls()
  out <- filter_allele_table(calls)
  expect_equal(out, calls)
})

test_that("filter_allele_table keeps extra columns", {
  calls <- dplyr::mutate(example_calls(), panel = "p1")
  out <- filter_allele_table(calls, min_reads = 20)
  expect_true("panel" %in% names(out))
  expect_false("wsaf" %in% names(out))
})

test_that("filter_allele_table validates columns, thresholds, and overwrite", {
  expect_error(
    filter_allele_table(data.frame(specimen_name = "s")),
    "Missing required columns"
  )
  expect_error(
    filter_allele_table(example_calls(), min_reads = -1),
    "min_reads"
  )
  expect_error(
    filter_allele_table(example_calls(), min_wsaf = 1.2),
    "min_wsaf"
  )
  expect_error(
    filter_allele_table(example_calls(), min_reads = 1000),
    "All alleles were removed"
  )

  tmp <- tempfile(fileext = ".tsv")
  on.exit(unlink(tmp), add = TRUE)
  readr::write_tsv(example_calls(), tmp)
  expect_error(
    filter_allele_table(example_calls(), output = tmp, overwrite = FALSE),
    "already exists"
  )
})

test_that("filter_allele_table works on packaged example", {
  path <- system.file(
    "extdata",
    "example_allele_table.tsv",
    package = "PGEcore"
  )
  skip_if_not(nzchar(path) && file.exists(path))
  out <- filter_allele_table(path, min_reads = 50)
  expect_true(nrow(out) >= 1)
  expect_true(all(out$reads >= 50))
  expect_true(all(
    c("specimen_name", "target_name", "seq", "reads") %in% names(out)
  ))
})
