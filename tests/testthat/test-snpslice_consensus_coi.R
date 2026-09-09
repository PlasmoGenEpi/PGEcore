# Chains shaped the way snpslice_consensus_coi() reads them: each carries the
# MAP allocation (hosts x strains) and dictionary (strains x loci) matrices.
make_chain <- function(alloc, dict) {
  list(map_allocation_matrix = alloc, map_dictionary_matrix = dict)
}

test_that("a haplotype carried by every restart keeps close to its raw count", {
  dict <- rbind(c(1, 0), c(0, 1))
  alloc <- rbind(c(1, 0), c(1, 1))
  chains <- list(make_chain(alloc, dict), make_chain(alloc, dict))
  out <- snpslice_consensus_coi(chains, "map")
  expect_length(out, 2L)
  expect_true(all(out >= 1))
  expect_gt(out[2], out[1])
})

test_that("a haplotype seen in one restart of four is discounted", {
  dict <- rbind(c(1, 0), c(0, 1))
  shared <- rbind(c(1, 0), c(1, 0))
  extra <- rbind(c(1, 1), c(1, 0))
  steady <- list(make_chain(shared, dict), make_chain(shared, dict),
                 make_chain(shared, dict), make_chain(shared, dict))
  sporadic <- list(make_chain(extra, dict), make_chain(shared, dict),
                   make_chain(shared, dict), make_chain(shared, dict))
  expect_lt(
    snpslice_consensus_coi(sporadic, "map")[1],
    sum(extra[1, ])
  )
  expect_equal(
    snpslice_consensus_coi(steady, "map")[2],
    snpslice_consensus_coi(sporadic, "map")[2]
  )
})

test_that("duplicate dictionary rows are collapsed before counting", {
  dict_dup <- rbind(c(1, 0), c(1, 0))
  dict_uniq <- rbind(c(1, 0), c(0, 1))
  dup <- snpslice_consensus_coi(
    list(make_chain(rbind(c(1, 1), c(1, 0)), dict_dup)), "map")
  uniq <- snpslice_consensus_coi(
    list(make_chain(rbind(c(1, 1), c(1, 0)), dict_uniq)), "map")
  expect_lt(dup[1], uniq[1])
})

test_that("every host is floored at one strain", {
  dict <- rbind(c(1, 0), c(0, 1))
  alloc <- rbind(c(0, 0), c(1, 1))
  out <- snpslice_consensus_coi(list(make_chain(alloc, dict)), "map")
  expect_equal(out[1], 1)
})

test_that("a bare single-chain result is accepted", {
  dict <- rbind(c(1, 0), c(0, 1))
  alloc <- rbind(c(1, 0), c(1, 1))
  expect_equal(
    snpslice_consensus_coi(make_chain(alloc, dict), "map"),
    snpslice_consensus_coi(list(make_chain(alloc, dict)), "map")
  )
})
