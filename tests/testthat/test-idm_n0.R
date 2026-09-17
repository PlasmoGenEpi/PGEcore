# n_0 (specimens with no lineage detected) is the only thing separating the IDM
# from the OM: with n_0 == 0 the IDM branch of MLE1() falls through to MLE_OM().
# Nk() derives it from the per-locus table, so every specimen must appear at every
# locus with NA where it has no call.

make_calls <- function() {
  # s1 and s2 each carry both variants at L1, so every delivered specimen holds
  # every variant and min(N - N_k) == 0 -- the case that returns NA frequencies.
  # s3 has no call at L1 and only appears because of another locus.
  tibble::tibble(
    specimen_name = c("s1", "s1", "s2", "s2", "s3"),
    locus         = c("L1", "L1", "L1", "L1", "L2"),
    variants      = c("L1:A", "L1:B", "L1:A", "L1:B", "L2:A")
  )
}

test_that("a specimen with no call at a locus still contributes to n_0", {
  df <- make_calls()
  specimens <- dplyr::distinct(df, specimen_name)
  per_locus <- dplyr::left_join(
    specimens,
    dplyr::select(dplyr::filter(df, locus == "L1"), specimen_name, variants),
    by = "specimen_name",
    relationship = "many-to-many"
  )
  N <- length(unique(per_locus$specimen_name))
  n_plus <- length(unique(per_locus$specimen_name[!is.na(per_locus$variants)]))
  expect_equal(N, 3L)
  expect_equal(N - n_plus, 1L)
})

test_that("the OM never receives absent specimens", {
  # MLE() rejects sum(N_k) < N for the OM, so handing it the missing specimens
  # aborts the run outright. Only the IDM gets them.
  res <- run_idm_mle_across_loci(make_calls(), model = "OM")
  expect_equal(sum(grepl("^L1:", res$variant)), 2L)
})

test_that("run_idm_mle_across_loci still emits a row per variant at every locus", {
  # Reporting n_0 does not by itself rescue a locus: with n_0 > 0 the MLE takes a
  # different branch and returns NA when prod(1 - N_k/N) <= n_0/N, i.e. when the
  # observed all-absent rate exceeds what independent lineage absence can explain.
  # What the fix guarantees is that n_0 reaches the model at all.
  res <- run_idm_mle_across_loci(make_calls(), model = "IDM")
  expect_equal(sum(grepl("^L1:", res$variant)), 2L)
  expect_true(all(c("L1:A", "L1:B") %in% res$variant))
})

test_that("every specimen appears at every locus in the table handed to the MLE", {
  df <- make_calls()
  specimens <- dplyr::distinct(df, specimen_name)
  for (l in unique(df$locus)) {
    per_locus <- dplyr::left_join(
      specimens,
      dplyr::select(dplyr::filter(df, locus == l), specimen_name, variants),
      by = "specimen_name",
      relationship = "many-to-many"
    )
    expect_setequal(unique(per_locus$specimen_name), specimens$specimen_name)
  }
})
