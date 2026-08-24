test_that("estimate_coi_naive_from_alleles integer method picks nth count", {
  alleles <- tibble::tibble(
    specimen_name = c("s1", "s1", "s1", "s1", "s1"),
    target_name = c("L1", "L1", "L2", "L3", "L3"),
    reads = c(10L, 20L, 15L, 5L, 8L),
    seq = c("A", "B", "A", "A", "C")
  )
  # L1: 2 alleles, L2: 1, L3: 2 -> sorted desc: 2, 2, 1
  out <- estimate_coi_naive_from_alleles(
    alleles,
    method = "integer_method",
    integer_threshold = 1
  )
  expect_equal(out$specimen_name, "s1")
  expect_equal(out$coi, 2)

  out2 <- estimate_coi_naive_from_alleles(
    alleles,
    method = "integer_method",
    integer_threshold = 3
  )
  expect_equal(out2$coi, 1)
})

test_that("estimate_coi_naive_from_alleles quantile method uses n_limit", {
  alleles <- tibble::tibble(
    specimen_name = rep("s1", 4),
    target_name = c("L1", "L2", "L3", "L4"),
    reads = c(10L, 20L, 15L, 5L),
    seq = c("A", "B", "C", "D")
  )
  # one allele per locus; loci (legacy row count) = 4
  # n_limit = floor((4-1)*0) + 1 = 1 -> first (highest) = 1
  out <- estimate_coi_naive_from_alleles(
    alleles,
    method = "quantile_method",
    quantile_threshold = 0
  )
  expect_equal(out$coi, 1)
})

test_that("estimate_coi_naive validates method and thresholds", {
  alleles <- tibble::tibble(
    specimen_name = "s1",
    target_name = "L1",
    reads = 10L,
    seq = "A"
  )
  expect_error(
    estimate_coi_naive_from_alleles(alleles, method = "nope"),
    regexp = "."
  )
  expect_error(
    estimate_coi_naive_from_alleles(
      alleles,
      method = "integer_method",
      integer_threshold = 0
    ),
    regexp = "."
  )
  expect_error(
    estimate_coi_naive_from_alleles(
      alleles,
      method = "quantile_method",
      quantile_threshold = 1.5
    ),
    regexp = "."
  )
})

test_that("estimate_coi_naive works on packaged example allele table", {
  path <- system.file(
    "extdata",
    "example_allele_table.tsv",
    package = "PGEcore"
  )
  skip_if_not(nzchar(path) && file.exists(path))

  out <- estimate_coi_naive(path, method = "integer_method", integer_threshold = 1)
  expect_true(all(c("specimen_name", "coi") %in% names(out)))
  expect_true(nrow(out) > 0)
  expect_true(all(out$coi >= 1))

  tmp_out <- tempfile(fileext = ".tsv")
  on.exit(unlink(tmp_out), add = TRUE)
  out2 <- estimate_coi_naive(
    path,
    output_path = tmp_out,
    method = "integer_method"
  )
  expect_true(file.exists(tmp_out))
  expect_equal(out2$coi, out$coi)
})
