# Fake MOIRe chain shaped the way extract_moire_chain_draws() reads it:
# per-sample draw vectors, per-locus lists of per-iteration allele-frequency
# vectors, and a mean_coi vector.
make_chain <- function(n_iter, n_samples = 2, n_loci = 1, n_alleles = 2) {
  per_sample <- function() {
    lapply(seq_len(n_samples), function(i) stats::runif(n_iter))
  }
  list(
    coi = per_sample(),
    eps_pos = per_sample(),
    eps_neg = per_sample(),
    relatedness = per_sample(),
    allele_freqs = lapply(seq_len(n_loci), function(l) {
      lapply(seq_len(n_iter), function(it) stats::runif(n_alleles))
    }),
    mean_coi = stats::runif(n_iter)
  )
}

truncate_chain <- function(chain, k) {
  keep <- function(v) v[seq_len(k)]
  chain$coi <- lapply(chain$coi, keep)
  chain$eps_pos <- lapply(chain$eps_pos, keep)
  chain$eps_neg <- lapply(chain$eps_neg, keep)
  chain$relatedness <- lapply(chain$relatedness, keep)
  chain$allele_freqs <- lapply(chain$allele_freqs, keep)
  chain$mean_coi <- keep(chain$mean_coi)
  chain
}

fake_results <- function(chains) {
  list(
    args = list(data = list(sample_ids = c("s1", "s2"), loci = "L1")),
    chains = chains
  )
}

test_that("prepare_moire_convergence_output summarizes equal-length chains", {
  skip_if_not_installed("posterior")
  set.seed(1)
  out <- prepare_moire_convergence_output(
    fake_results(list(make_chain(120), make_chain(120)))
  )
  # 2 samples x (coi, eps_pos, eps_neg, relatedness) + 2 alleles + mean_coi
  expect_equal(nrow(out), 11)
  expect_true(all(c("variable", "rhat", "ess_bulk") %in% names(out)))
  expect_true(all(is.finite(out$rhat)))
})

test_that("prepare_moire_convergence_output truncates ragged chains to the shortest", {
  skip_if_not_installed("posterior")
  set.seed(2)
  full <- list(make_chain(120), make_chain(120))
  ragged <- list(full[[1]], truncate_chain(full[[2]], 40))

  expect_warning(
    got <- prepare_moire_convergence_output(fake_results(ragged)),
    "unequal numbers of draws"
  )
  # Truncating every chain up front must give exactly the same diagnostics,
  # i.e. the fix keeps each chain's FIRST n_iter draws rather than its tail.
  want <- prepare_moire_convergence_output(
    fake_results(lapply(full, truncate_chain, 40))
  )
  expect_equal(got, want)
  expect_equal(nrow(got), 11)
})

test_that("prepare_moire_convergence_output errors when a chain recorded no draws", {
  skip_if_not_installed("posterior")
  set.seed(3)
  expect_error(
    prepare_moire_convergence_output(
      fake_results(list(make_chain(120), truncate_chain(make_chain(120), 0)))
    ),
    "recorded no draws for 1 of 2 chains"
  )
})

test_that("assert_moire_chains_have_draws passes intact chains and counts empty ones", {
  set.seed(4)
  expect_silent(
    assert_moire_chains_have_draws(fake_results(list(make_chain(10), make_chain(10))))
  )
  # A cap landing in burn-in empties every chain, not just one.
  expect_error(
    assert_moire_chains_have_draws(
      fake_results(lapply(list(make_chain(10), make_chain(10)), truncate_chain, 0))
    ),
    "recorded no draws for 2 of 2 chains"
  )
  # The message has to name the actual remedy, since the caller cannot tell
  # from moire's own error what went wrong.
  expect_error(
    assert_moire_chains_have_draws(
      fake_results(list(make_chain(10), truncate_chain(make_chain(10), 0)))
    ),
    "max_runtime"
  )
})
