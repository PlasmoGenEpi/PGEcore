# Fake snp_slice() result shaped the way prepare_snpslice_optim_output() reads
# it: a $chains list whose elements carry a $diagnostics record and a
# $map_allocation_matrix of hosts x strains.
make_chain <- function(chain_id, seed, map_logpost, alloc,
                       map_iteration = 80L, final_iteration = 100L,
                       map_kstar = 3L, map_ktrunc = 4L) {
  list(
    diagnostics = list(
      chain_id = chain_id,
      seed = seed,
      map_logpost = map_logpost,
      map_iteration = map_iteration,
      final_iteration = final_iteration,
      map_kstar = map_kstar,
      map_ktrunc = map_ktrunc
    ),
    map_allocation_matrix = alloc
  )
}

# Hosts x strains allocation with a given per-host strain count.
alloc_from_coi <- function(coi, n_strains = 5L) {
  t(vapply(coi, function(k) {
    row <- rep(0L, n_strains)
    row[seq_len(k)] <- 1L
    row
  }, integer(n_strains)))
}

test_that("prepare_snpslice_optim_output returns one row per restart", {
  coi <- c(1L, 2L, 3L, 4L, 5L)
  res <- list(
    chains = list(
      make_chain(1L, 11L, -100, alloc_from_coi(coi)),
      make_chain(2L, 22L, -90, alloc_from_coi(coi)),
      make_chain(3L, 33L, -110, alloc_from_coi(coi))
    ),
    best_chain = 2L
  )
  out <- prepare_snpslice_optim_output(res)

  expect_s3_class(out, "tbl_df")
  expect_equal(nrow(out), 3L)
  expect_equal(
    colnames(out),
    c("chain_id", "seed", "map_logpost", "is_best", "map_iteration",
      "final_iteration", "plateau_frac", "map_kstar", "map_ktrunc",
      "coi_mean", "coi_ccc_to_best")
  )
  expect_equal(out$chain_id, c(1L, 2L, 3L))
  expect_equal(out$seed, c(11L, 22L, 33L))
  expect_equal(out$map_logpost, c(-100, -90, -110))
})

test_that("is_best marks exactly the chain snp_slice() reported", {
  coi <- c(2L, 2L, 3L, 3L, 4L)
  res <- list(
    chains = list(
      make_chain(1L, 11L, -100, alloc_from_coi(coi)),
      make_chain(2L, 22L, -90, alloc_from_coi(coi)),
      make_chain(3L, 33L, -110, alloc_from_coi(coi))
    ),
    best_chain = 2L
  )
  out <- prepare_snpslice_optim_output(res)

  expect_equal(sum(out$is_best), 1L)
  expect_true(out$is_best[2])
  # best_chain is an index into $chains, not a chain_id lookup
  expect_equal(out$chain_id[out$is_best], 2L)
})

test_that("coi_ccc_to_best is 1 for the reported chain and lower for a differing one", {
  best_coi <- c(1L, 2L, 3L, 4L, 5L)
  other_coi <- c(5L, 4L, 3L, 2L, 1L)
  res <- list(
    chains = list(
      make_chain(1L, 11L, -90, alloc_from_coi(best_coi)),
      make_chain(2L, 22L, -100, alloc_from_coi(other_coi))
    ),
    best_chain = 1L
  )
  out <- prepare_snpslice_optim_output(res)

  expect_equal(out$coi_ccc_to_best[1], 1)
  expect_lt(out$coi_ccc_to_best[2], out$coi_ccc_to_best[1])
})

test_that("coi_mean is the mean per-host strain count", {
  res <- list(
    chains = list(make_chain(1L, 11L, -90, alloc_from_coi(c(1L, 2L, 3L)))),
    best_chain = 1L
  )
  out <- prepare_snpslice_optim_output(res)
  expect_equal(out$coi_mean, 2)
})

test_that("plateau_frac is map_iteration over final_iteration", {
  res <- list(
    chains = list(
      make_chain(1L, 11L, -90, alloc_from_coi(c(1L, 2L, 3L)),
                 map_iteration = 25L, final_iteration = 100L)
    ),
    best_chain = 1L
  )
  out <- prepare_snpslice_optim_output(res)
  expect_equal(out$plateau_frac, 0.25)
})

test_that("a single-chain result with no $chains is handled", {
  single <- make_chain(1L, 11L, -90, alloc_from_coi(c(2L, 2L, 4L)))
  out <- prepare_snpslice_optim_output(single)

  expect_equal(nrow(out), 1L)
  expect_true(out$is_best)
  expect_equal(out$coi_ccc_to_best, 1)
  expect_equal(out$coi_mean, 8 / 3)
})

test_that("missing chain_id and seed fall back to the position and NA", {
  chain <- make_chain(1L, 11L, -90, alloc_from_coi(c(1L, 2L, 3L)))
  chain$diagnostics$chain_id <- NULL
  chain$diagnostics$seed <- NULL
  res <- list(chains = list(chain), best_chain = 1L)
  out <- prepare_snpslice_optim_output(res)

  expect_equal(out$chain_id, 1L)
  expect_true(is.na(out$seed))
})

test_that("snpslice_ccc handles degenerate input", {
  expect_true(is.na(snpslice_ccc(c(1, 2), c(1, 2))))
  expect_equal(snpslice_ccc(rep(3, 5), rep(3, 5)), 1)
  expect_equal(snpslice_ccc(c(1, 2, 3, 4), c(1, 2, 3, 4)), 1)
  expect_true(is.na(snpslice_ccc(c(1, NA, 3), c(1, 2, NA))))
})
