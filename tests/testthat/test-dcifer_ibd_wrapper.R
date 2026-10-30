test_that("dcifer_pair_chunks covers every pair exactly once", {
  for (pairs in 1:7) {
    for (cores in 1:5) {
      idx <- unlist(dcifer_pair_chunks(pairs, cores))
      expect_identical(sort(idx), seq_len(pairs), info = paste(pairs, cores))
      expect_true(all(lengths(dcifer_pair_chunks(pairs, cores)) > 0))
    }
  }
})

test_that("dcifer_test_side follows the relatedness null", {
  expect_equal(dcifer_test_side(0), "right")
  expect_equal(dcifer_test_side(0L), "right")
  expect_equal(dcifer_test_side(1), "left")
  expect_equal(dcifer_test_side(0.25), "two-sided")
})

test_that("dcifer_ibd_wrapper rejects allele_freq_table with pop_name_col", {
  skip_if_not_installed("dcifer")
  skip_if_not_installed("foreach")
  skip_if_not_installed("doParallel")
  skip_if_not_installed("parallelly")
  skip_if_not_installed("iterators")
  expect_error(
    dcifer_ibd_wrapper(
      allele_table = "a.tsv", relatedness_output = "o.tsv",
      allele_freq_table = "f.tsv", pop_name_col = "population"
    ),
    "cannot be used with --pop_name_col"
  )
})

# Two specimens, listed with "zed" before "abe" so the table is not in the
# alphabetical order dcifer uses internally; specimen names carry dots.
make_dcifer_table <- function() {
  data.frame(
    specimen_name = c("zed.1", "zed.1", "zed.1", "abe.1", "abe.1", "kim.1", "kim.1"),
    target_name   = c("L1", "L1", "L2", "L1", "L2", "L1", "L2"),
    seq           = c("a", "b", "x", "a", "x", "b", "y"),
    reads         = 10L
  )
}

test_that("supplied allele frequencies are aligned to dcifer's allele order", {
  skip_if_not_installed("dcifer")
  tab <- make_dcifer_table()
  dsmp <- dcifer::formatDat(
    dplyr::rename(tab, target_value = "seq"),
    svar = "specimen_name", lvar = "target_name", avar = "target_value"
  )
  # frequencies listed with loci and alleles in reverse order
  freq_path <- tempfile(fileext = ".tsv")
  on.exit(unlink(freq_path), add = TRUE)
  readr::write_tsv(data.frame(
    target_name = c("L2", "L2", "L1", "L1"),
    seq         = c("y", "x", "b", "a"),
    freq        = c(0.1, 0.9, 0.3, 0.7)
  ), freq_path)
  afreq <- create_dcifer_allele_freq_input(freq_path, dsmp)
  expect_identical(names(afreq), names(dsmp[[1]]))
  expect_identical(names(afreq$L1), names(dsmp[[1]]$L1))
  expect_identical(names(afreq$L2), names(dsmp[[1]]$L2))
  expect_equal(unname(afreq$L1[c("a", "b")]), c(0.7, 0.3))
  expect_equal(unname(afreq$L2[c("x", "y")]), c(0.9, 0.1))
})

test_that("a single population gives the same relatedness with or without pop_name_col", {
  skip_if_not_installed("dcifer")
  skip_if_not_installed("foreach")
  skip_if_not_installed("doParallel")
  skip_if_not_installed("parallelly")
  skip_if_not_installed("iterators")
  # "zed" is polyclonal and listed first, so the table is not in the
  # alphabetical order dcifer uses internally
  tab <- data.frame(
    specimen_name = c(rep("zed.1", 8), rep("abe.1", 4), rep("kim.1", 4), rep("joe.1", 4)),
    target_name   = c(rep(c("L1", "L2", "L3", "L4"), each = 2), rep(c("L1", "L2", "L3", "L4"), 3)),
    seq           = c("a", "b", "x", "y", "p", "q", "m", "n",
                      "a", "x", "p", "m",
                      "b", "y", "q", "n",
                      "a", "y", "p", "n"),
    reads         = 10L
  )
  coi <- data.frame(specimen_name = c("abe.1", "joe.1", "kim.1", "zed.1"), coi = c(1L, 1L, 1L, 3L))
  meta <- data.frame(specimen_name = coi$specimen_name, population = "p1")
  d <- tempfile(); dir.create(d)
  on.exit(unlink(d, recursive = TRUE), add = TRUE)
  readr::write_tsv(tab, file.path(d, "a.tsv"))
  readr::write_tsv(coi, file.path(d, "coi.tsv"))
  readr::write_tsv(meta, file.path(d, "meta.tsv"))
  run <- function(...) {
    res <- dcifer_ibd_wrapper(
      allele_table = file.path(d, "a.tsv"),
      relatedness_output = file.path(d, "o.tsv"),
      coi_table = file.path(d, "coi.tsv"),
      threads = 1L,
      ...
    )
    res[order(res$sample_a, res$sample_b), c("sample_a", "sample_b", "estimate")]
  }
  plain <- run()
  by_pop <- run(specimen_metadata = file.path(d, "meta.tsv"), pop_name_col = "population")
  expect_equal(nrow(plain), 6)
  expect_equal(by_pop, plain, ignore_attr = TRUE)
})

test_that("dcifer_ibd_wrapper infers COI for dotted specimen names", {
  skip_if_not_installed("dcifer")
  skip_if_not_installed("foreach")
  skip_if_not_installed("doParallel")
  skip_if_not_installed("parallelly")
  skip_if_not_installed("iterators")
  d <- tempfile(); dir.create(d)
  on.exit(unlink(d, recursive = TRUE), add = TRUE)
  readr::write_tsv(make_dcifer_table(), file.path(d, "a.tsv"))
  res <- dcifer_ibd_wrapper(
    allele_table = file.path(d, "a.tsv"),
    relatedness_output = file.path(d, "o.tsv"),
    threads = 1L
  )
  expect_equal(nrow(res), 3)
  expect_false(any(duplicated(res[, c("sample_a", "sample_b")])))
})

test_that("a single pair over several workers is reported once", {
  skip_if_not_installed("dcifer")
  skip_if_not_installed("foreach")
  skip_if_not_installed("doParallel")
  skip_if_not_installed("parallelly")
  skip_if_not_installed("iterators")
  tab <- make_dcifer_table()
  tab <- tab[tab$specimen_name != "kim.1", ]
  d <- tempfile(); dir.create(d)
  on.exit(unlink(d, recursive = TRUE), add = TRUE)
  readr::write_tsv(tab, file.path(d, "a.tsv"))
  res <- dcifer_ibd_wrapper(
    allele_table = file.path(d, "a.tsv"),
    relatedness_output = file.path(d, "o.tsv"),
    threads = 2L
  )
  expect_equal(nrow(res), 1)
})

test_that("estm_mmax and estm_equalr are passed to ibdEstM", {
  skip_if_not_installed("dcifer")
  skip_if_not_installed("foreach")
  skip_if_not_installed("doParallel")
  skip_if_not_installed("parallelly")
  skip_if_not_installed("iterators")
  skip_if_not_installed("checkmate")
  d <- tempfile(); dir.create(d)
  on.exit(unlink(d, recursive = TRUE), add = TRUE)
  readr::write_tsv(make_dcifer_table(), file.path(d, "a.tsv"))
  run <- function(...) {
    dcifer_ibd_wrapper(
      allele_table = file.path(d, "a.tsv"),
      relatedness_output = file.path(d, "o.tsv"),
      use_estm = TRUE,
      threads = 1L,
      ...
    )
  }
  expect_lte(max(run(estm_mmax = 1)$strain_pair), 1)
  expect_no_error(run(estm_equalr = TRUE))
  expect_error(run(estm_mmax = 7), "not supported")
})
