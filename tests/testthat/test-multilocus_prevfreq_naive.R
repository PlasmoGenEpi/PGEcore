test_that("calculate_multilocus_af_prev_presence_absence uses counts", {
  calls <- tibble::tibble(
    specimen_name = c("s1", "s1", "s2"),
    variant = c("A", "B", "A")
  )
  out <- calculate_multilocus_af_prev_presence_absence(calls)
  a <- out[out$variant == "A", ]
  b <- out[out$variant == "B", ]
  expect_equal(a$prev, 1)
  expect_equal(a$freq, 2 / 3)
  expect_equal(b$prev, 0.5)
  expect_equal(b$freq, 1 / 3)
  expect_equal(a$sample_total, 2)
  expect_equal(a$allele_total, 3)
})

test_that("calculate_multilocus_af_prev_wsaf_prop weights by wsaf", {
  calls <- tibble::tibble(
    specimen_name = c("s1", "s1", "s2"),
    variant = c("A", "B", "A"),
    wsaf = c(0.7, 0.3, 1)
  )
  out <- calculate_multilocus_af_prev_wsaf_prop(calls)
  a <- out[out$variant == "A", ]
  b <- out[out$variant == "B", ]
  expect_equal(a$prev, 1)
  expect_equal(a$freq, 1.7 / 2)
  expect_equal(b$prev, 0.5)
  expect_equal(b$freq, 0.3 / 2)
})

test_that("generate_single_locus_prev_freq_from_multilocus_groups_wsaf_prop splits variants", {
  calls <- tibble::tibble(
    group_id = "g1",
    specimen_name = "s1",
    variant = "PF3D7_0417200.1:51:N;PF3D7_0810800.1:437:A",
    wsaf = 1
  )
  out <- generate_single_locus_prev_freq_from_multilocus_groups_wsaf_prop(calls)
  expect_equal(nrow(out), 2)
  expect_true(all(out$freq == 1))
  expect_true(all(out$prev == 1))
  expect_true(all(grepl(":", out$variant)))
})

test_that("multilocus_prevfreq_naive rejects invalid method", {
  aa_path <- system.file(
    "extdata",
    "example2_amino_acid_calls.tsv",
    package = "PGEcore"
  )
  groups_path <- system.file(
    "extdata",
    "example_loci_groups.tsv",
    package = "PGEcore"
  )
  skip_if_not(nzchar(aa_path) && file.exists(aa_path))
  skip_if_not(nzchar(groups_path) && file.exists(groups_path))

  expect_error(
    multilocus_prevfreq_naive(aa_path, groups_path, method = "nope"),
    "is not a valid method"
  )
})

test_that("multilocus_prevfreq_naive works on packaged example data", {
  aa_path <- system.file(
    "extdata",
    "example2_amino_acid_calls.tsv",
    package = "PGEcore"
  )
  groups_path <- system.file(
    "extdata",
    "example_loci_groups.tsv",
    package = "PGEcore"
  )
  skip_if_not(nzchar(aa_path) && file.exists(aa_path))
  skip_if_not(nzchar(groups_path) && file.exists(groups_path))

  out <- multilocus_prevfreq_naive(aa_path, groups_path)
  expect_true(nrow(out) > 0)
  expect_true(all(c("variant", "prev", "freq", "group_id") %in% names(out)))
  expect_true(all(out$prev >= 0 & out$prev <= 1))
  expect_true(all(out$freq >= 0 & out$freq <= 1))

  pa <- multilocus_prevfreq_naive(
    aa_path,
    groups_path,
    method = "presence_absence"
  )
  expect_true(nrow(pa) > 0)
  expect_true("allele_count" %in% names(pa))
  expect_false("allele_count" %in% names(out))

  tmp_out <- tempfile(fileext = ".tsv")
  tmp_sl <- tempfile(fileext = ".tsv")
  on.exit(unlink(c(tmp_out, tmp_sl)), add = TRUE)
  written <- multilocus_prevfreq_naive(
    aa_path,
    groups_path,
    output_path = tmp_out,
    recalc_single_locus_output_path = tmp_sl
  )
  expect_true(file.exists(tmp_out))
  expect_true(file.exists(tmp_sl))
  expect_equal(nrow(written), nrow(out))
})
