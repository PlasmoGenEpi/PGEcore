# The OM's score equation is
#   f(lambda) = lambda + sum_k log(1 - n_k (1 - exp(-lambda)))
# A lineage at prevalence exactly 1 contributes log(exp(-lambda)) = -lambda,
# cancelling the leading term and leaving f strictly negative: no positive root,
# so the MLE does not exist and the published code returns NA. "OM_CC" applies a
# continuity correction to n_k, which holds it off that boundary.

# 40 specimens; L1:A is carried by every one of them.
fixed_variant_calls <- function() {
  dplyr::bind_rows(
    tibble::tibble(
      specimen_name = rep(sprintf("s%02d", 1:40), each = 2),
      locus = "L1",
      variants = rep(c("L1:A", "L1:B"), times = 40)
    ) |> dplyr::filter(!(specimen_name %in% sprintf("s%02d", 31:40) &
                           variants == "L1:B")),
    tibble::tibble(specimen_name = "s01", locus = "L2", variants = "L2:A")
  )
}

l1 <- function(res) res$freq[grepl("^L1:", res$variant)]

test_that("the OM cannot solve a locus with a lineage at prevalence 1", {
  res <- suppressWarnings(
    run_idm_mle_across_loci(fixed_variant_calls(), model = "OM")
  )
  expect_true(all(is.na(l1(res))))
})

test_that("OM_CC solves it, and ranks the near-fixed lineage highest", {
  res <- suppressWarnings(
    run_idm_mle_across_loci(fixed_variant_calls(), model = "OM_CC")
  )
  freqs <- l1(res)
  expect_false(any(is.na(freqs)))
  expect_true(all(freqs >= 0 & freqs <= 1))
  variants <- res$variant[grepl("^L1:", res$variant)]
  expect_gt(freqs[variants == "L1:A"], freqs[variants == "L1:B"])
})

# No lineage at prevalence 1 here, so the published estimator has a root.
solvable_calls <- function() {
  dplyr::bind_rows(
    tibble::tibble(
      specimen_name = rep(sprintf("s%02d", 1:30), each = 2),
      locus = "L1",
      variants = rep(c("L1:A", "L1:B"), times = 30)
    ) |> dplyr::filter(!(specimen_name %in% sprintf("s%02d", 1:12) &
                           variants == "L1:B")),
    tibble::tibble(
      specimen_name = sprintf("s%02d", 31:40), locus = "L1", variants = "L1:B"
    ),
    tibble::tibble(specimen_name = "s01", locus = "L2", variants = "L2:A")
  )
}

test_that("OM_CC leaves loci the OM already solves essentially untouched", {
  df <- solvable_calls()
  om <- suppressWarnings(run_idm_mle_across_loci(df, model = "OM"))
  cc <- suppressWarnings(run_idm_mle_across_loci(df, model = "OM_CC"))
  expect_false(any(is.na(l1(om))))
  expect_false(any(is.na(l1(cc))))
  expect_lt(max(abs(l1(om) - l1(cc))), 0.05)
})

test_that("the continuity knob defaults to the published estimator", {
  # MLE(continuity = 0) must be the untouched published path, so "OM" and an
  # explicit continuity of 0 have to agree exactly -- checked on a locus that
  # actually solves, so this compares numbers rather than two NAs.
  df <- solvable_calls()
  plain <- suppressWarnings(run_idm_mle_across_loci(df, model = "OM"))
  explicit <- suppressWarnings(
    idm_locus_mle(
      dplyr::select(dplyr::filter(df, locus == "L1"), specimen_name, variants),
      "OM", 1.0, 0.1, continuity = 0
    )
  )
  expect_false(any(is.na(explicit$freq)))
  # The vendor hands back a 1-row matrix; the caller flattens it.
  expect_equal(l1(plain), as.numeric(explicit$freq))
})

test_that("IDM_wrapper accepts OM_CC", {
  expect_error(
    IDM_wrapper(aa_calls = "x", slaf_output = "y", model = "CC_OM"),
    "IDM \\| OM \\| IDM_OM \\| OM_CC"
  )
})

test_that("IDM_OM_CC falls back to the corrected OM, not the plain one", {
  # The IDM cannot solve this locus, and neither can the plain OM -- a lineage
  # sits at prevalence 1. Only the corrected fallback returns frequencies.
  df <- fixed_variant_calls()
  idm_om <- suppressWarnings(run_idm_mle_across_loci(df, model = "IDM_OM"))
  idm_om_cc <- suppressWarnings(run_idm_mle_across_loci(df, model = "IDM_OM_CC"))
  om_cc <- suppressWarnings(run_idm_mle_across_loci(df, model = "OM_CC"))

  expect_true(all(is.na(l1(idm_om))))
  expect_false(any(is.na(l1(idm_om_cc))))
  expect_equal(l1(idm_om_cc), l1(om_cc))
})

test_that("IDM_OM_CC leaves the IDM's own solves uncorrected", {
  # Where the IDM solves a locus the fallback must not fire, so the correction
  # never touches it and the answer matches plain IDM exactly.
  df <- solvable_calls()
  idm <- suppressWarnings(run_idm_mle_across_loci(df, model = "IDM"))
  hybrid <- suppressWarnings(run_idm_mle_across_loci(df, model = "IDM_OM_CC"))
  expect_false(any(is.na(l1(idm))))
  expect_equal(l1(hybrid), l1(idm))
})

test_that("IDM_wrapper accepts IDM_OM_CC", {
  expect_error(
    IDM_wrapper(aa_calls = "x", slaf_output = "y", model = "CC_IDM_OM"),
    "IDM \\| OM \\| IDM_OM \\| OM_CC \\| IDM_OM_CC"
  )
})
