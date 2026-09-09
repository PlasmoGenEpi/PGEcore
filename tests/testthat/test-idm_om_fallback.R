# model = "IDM_OM" runs the IDM and re-solves with the OM at any locus the IDM
# leaves unsolved. The IDM returns NA when prod(1 - N_k/N) <= n_0/N, which a
# near-fixed variant reaches easily: the product collapses toward zero while any
# missing specimen keeps n_0/N above it.

# 40 specimens, both variants near-fixed (37 carry both), one specimen absent.
# prod(1 - 38/40)^2 = 0.0025 <= n_0/N = 0.025, so the IDM returns NA. The OM sees
# only the 39 delivered specimens, where min(N - N_k) == 1, and solves.
near_fixed_calls <- function() {
  both <- sprintf("s%02d", 1:37)
  dplyr::bind_rows(
    tibble::tibble(
      specimen_name = rep(both, each = 2),
      locus = "L1",
      variants = rep(c("L1:A", "L1:B"), times = 37)
    ),
    tibble::tibble(
      specimen_name = c("s38", "s39"),
      locus = "L1",
      variants = c("L1:A", "L1:B")
    ),
    tibble::tibble(specimen_name = "s40", locus = "L2", variants = "L2:A")
  )
}

# Asymmetric variant counts far from fixation, so the IDM solves the locus and
# lands somewhere the OM does not.
solvable_calls <- function() {
  dplyr::bind_rows(
    tibble::tibble(
      specimen_name = sprintf("s%02d", 1:8), locus = "L1", variants = "L1:A"
    ),
    tibble::tibble(
      specimen_name = sprintf("s%02d", 9:11), locus = "L1", variants = "L1:B"
    ),
    tibble::tibble(
      specimen_name = rep(sprintf("s%02d", 12:18), each = 2),
      locus = "L1",
      variants = rep(c("L1:A", "L1:B"), times = 7)
    ),
    tibble::tibble(
      specimen_name = c("s19", "s20"), locus = "L2", variants = "L2:A"
    )
  )
}

l1 <- function(res) res$freq[grepl("^L1:", res$variant)]

test_that("IDM_OM recovers a locus the IDM leaves unsolved", {
  df <- near_fixed_calls()
  idm <- suppressWarnings(run_idm_mle_across_loci(df, model = "IDM"))
  om <- suppressWarnings(run_idm_mle_across_loci(df, model = "OM"))
  hybrid <- suppressWarnings(run_idm_mle_across_loci(df, model = "IDM_OM"))

  expect_true(all(is.na(l1(idm))))
  expect_false(any(is.na(l1(om))))
  expect_equal(l1(hybrid), l1(om))
})

test_that("IDM_OM leaves a locus the IDM does solve alone", {
  df <- solvable_calls()
  idm <- suppressWarnings(run_idm_mle_across_loci(df, model = "IDM"))
  om <- suppressWarnings(run_idm_mle_across_loci(df, model = "OM"))
  hybrid <- suppressWarnings(run_idm_mle_across_loci(df, model = "IDM_OM"))

  expect_false(any(is.na(l1(idm))))
  # The two models disagree here, so matching the IDM is evidence the fallback
  # did not fire rather than a coincidence of equal answers.
  expect_false(isTRUE(all.equal(l1(idm), l1(om))))
  expect_equal(l1(hybrid), l1(idm))
})

test_that("IDM_OM does not invent an answer when the OM also fails", {
  # Every delivered specimen carries every variant, so min(N - N_k) == 0 and both
  # models return NA. The fallback must pass the NA through.
  df <- tibble::tibble(
    specimen_name = c("s1", "s1", "s2", "s2", "s3"),
    locus = c("L1", "L1", "L1", "L1", "L2"),
    variants = c("L1:A", "L1:B", "L1:A", "L1:B", "L2:A")
  )
  hybrid <- suppressWarnings(run_idm_mle_across_loci(df, model = "IDM_OM"))
  expect_true(all(is.na(l1(hybrid))))
})

test_that("IDM_OM emits one row per variant at every locus", {
  hybrid <- suppressWarnings(
    run_idm_mle_across_loci(near_fixed_calls(), model = "IDM_OM")
  )
  expect_setequal(hybrid$variant, c("L1:A", "L1:B", "L2:A"))
})

test_that("IDM_wrapper accepts IDM_OM and rejects unknown models", {
  expect_error(
    IDM_wrapper(aa_calls = "x", slaf_output = "y", model = "OM_IDM"),
    "IDM \\| OM \\| IDM_OM"
  )
})
