# check_fem_biallelic() must drop a target with more than one alternate amino
# acid: FEM codes every alternate as 1, so two alternates would be conflated.
make_fem_calls <- function() {
  tibble::tibble(
    specimen_name = c("s1", "s2", "s3", "s1", "s2", "s1", "s2"),
    gene_id       = c("g", "g", "g", "g", "g", "g", "g"),
    aa_position   = c(1L, 1L, 1L, 2L, 2L, 3L, 3L),
    ref_aa        = c("A", "A", "A", "C", "C", "N", "N"),
    aa            = c("A", "K", "R", "C", "Y", "N", "N")
  ) |>
    dplyr::mutate(unique_targets = paste(gene_id, aa_position, sep = ":"))
}

test_that("check_fem_biallelic drops targets with two alternate amino acids", {
  calls <- make_fem_calls()
  expect_warning(res <- check_fem_biallelic(calls), "g:1")
  only_biallelic <- res[[1]]
  alt_calls <- res[[2]]
  monos <- res[[3]]
  # g:1 has alternates K and R: gone from both the kept calls and the alternates
  expect_false("g:1" %in% only_biallelic$unique_targets)
  expect_false("g:1" %in% alt_calls$unique_targets)
  # g:2 is a proper biallelic target and keeps its single alternate
  expect_true("g:2" %in% only_biallelic$unique_targets)
  expect_equal(unique(alt_calls$aa[alt_calls$unique_targets == "g:2"]), "Y")
  # g:3 is monoallelic and reported as such
  expect_equal(monos$unique_targets, "g:3")
  expect_false("g:3" %in% only_biallelic$unique_targets)
})

test_that("check_fem_biallelic is silent when every target is biallelic", {
  calls <- make_fem_calls() |> dplyr::filter(unique_targets == "g:2")
  expect_no_warning(res <- check_fem_biallelic(calls))
  expect_equal(unique(res[[1]]$unique_targets), "g:2")
})

# fem_seed_args() decides what reaches mcmc_sampling_parallel().
test_that("fem_seed_args forwards an explicit seed as an integer", {
  local_mocked_bindings(fem_supports_seed = function() TRUE)
  expect_identical(fem_seed_args(7), list(seed = 7L))
})

test_that("fem_seed_args draws a random seed when none is given", {
  local_mocked_bindings(fem_supports_seed = function() TRUE)
  a <- fem_seed_args(NULL)
  expect_true(is.integer(a$seed) && length(a$seed) == 1L && a$seed > 0L)
  b <- fem_seed_args(NULL)
  expect_false(identical(a$seed, b$seed))
})

test_that("fem_seed_args forwards nothing when the build has no seed argument", {
  local_mocked_bindings(fem_supports_seed = function() FALSE)
  expect_identical(fem_seed_args(NULL), list())
})

test_that("fem_seed_args errors on an explicit seed the build cannot take", {
  skip_if_not_installed("FreqEstimationModel")
  local_mocked_bindings(fem_supports_seed = function() FALSE)
  expect_error(fem_seed_args(1L), "no `seed` argument")
})

test_that("fem_supports_seed reads mcmc_sampling_parallel's formals", {
  skip_if_not_installed("FreqEstimationModel")
  expect_identical(
    fem_supports_seed(),
    "seed" %in% names(formals(FreqEstimationModel::mcmc_sampling_parallel))
  )
})
