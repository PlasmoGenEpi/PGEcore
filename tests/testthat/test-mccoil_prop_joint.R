# Read counts for 8 monoclonal, 8 two-strain (0.7/0.3) and 8 three-strain
# (0.6/0.28/0.12) samples at 40 loci, with strain proportions shared across
# loci and minor alleles below 10 reads zeroed as the pipeline does.
make_joint_snp_calls <- function(seed = 5) {
  set.seed(seed)
  k <- 40
  coi <- rep(1:3, each = 8)
  p <- stats::runif(k, 0.2, 0.8)
  props <- list(1, c(0.7, 0.3), c(0.6, 0.28, 0.12))
  rows <- list()
  for (i in seq_along(coi)) {
    w <- props[[coi[i]]]
    for (j in seq_len(k)) {
      s <- sum(w * stats::rbinom(length(w), 1, p[j]))
      a <- stats::rbinom(1, 300, s)
      b <- 300 - a
      if (min(a, b) < 10) {
        if (a < b) a <- 0 else b <- 0
      }
      rows[[length(rows) + 1]] <- data.frame(
        specimen_name = paste0("s", i), snp_name = paste0("l", j),
        seq_base = c("A", "T"), reads = c(a, b)
      )
    }
  }
  df <- do.call(rbind, rows)
  list(calls = df[df$reads > 0, ], coi = stats::setNames(coi, paste0("s", seq_along(coi))))
}

run_joint_wrapper <- function(sim, seed = 321) {
  f <- tempfile(fileext = ".tsv")
  readr::write_tsv(sim$calls, f)
  d <- tempfile()
  dir.create(d)
  THEREALMcCOIL_wrapper(
    f, file.path(d, "slaf.tsv"), file.path(d, "coi.tsv"),
    model = "proportional_joint", totalrun = 1500, burnin = 500,
    n_chains = 2, seed = seed,
    convergence_output = file.path(d, "conv.tsv"),
    convergence_summary_output = file.path(d, "convsum.tsv")
  )
}

test_that("proportional_joint rejects err_method other than 1", {
  df <- tibble::tibble(specimen_name = "S1", snp_name = "L1", seq_base = "A", reads = 5L)
  expect_error(
    PGEcore:::run_mccoil_chains(
      df, model = "proportional_joint", err_method = 3, work_dir = tempdir()
    ),
    "must be 1 for the proportional models"
  )
})

test_that("proportional_joint recovers COI on shared-proportion data", {
  skip_on_cran()
  sim <- make_joint_snp_calls()
  res <- run_joint_wrapper(sim)
  est <- stats::setNames(res$coi$coi, res$coi$specimen_name)[names(sim$coi)]
  expect_true(all(est[sim$coi == 1] == 1))
  expect_true(all(est[sim$coi == 2] == 2))
  expect_true(all(est[sim$coi == 3] >= 3 & est[sim$coi == 3] <= 4))
  # read noise and outlier rate are reported as extra parameters
  expect_true(all(c("rho", "pout") %in% res$convergence$variable))
  expect_equal(nrow(res$slaf), 40)
})

test_that("proportional_joint is reproducible for a given seed", {
  skip_on_cran()
  sim <- make_joint_snp_calls()
  a <- run_joint_wrapper(sim, seed = 7)
  b <- run_joint_wrapper(sim, seed = 7)
  expect_identical(a$coi, b$coi)
  expect_identical(a$slaf, b$slaf)
})

test_that("chains are identical whether run in parallel or one at a time", {
  skip_on_cran()
  skip_on_os("windows")
  df <- make_joint_snp_calls()$calls
  for (model in c("categorical", "proportional", "proportional_joint")) {
    run <- function(threads) {
      wd <- tempfile()
      PGEcore:::run_mccoil_chains(
        df, model = model, totalrun = 200, burnin = 50, n_chains = 3,
        threads = threads, seed = 11, work_dir = wd,
        threshold_ind = 5, threshold_site = 5
      )$traces
    }
    seq_tr <- run(1)
    par_tr <- run(3)
    for (i in 1:3) expect_identical(par_tr[[i]], seq_tr[[i]], info = model)
    # each chain has its own seed, so the chains themselves differ
    expect_false(identical(seq_tr[[1]], seq_tr[[2]]), info = model)
  }
})

test_that("proportional_joint chains start spread from 1 to 15 strains", {
  skip_on_cran()
  df <- make_joint_snp_calls()$calls
  tr <- PGEcore:::run_mccoil_chains(
    df, model = "proportional_joint", totalrun = 5, burnin = 1, n_chains = 3,
    seed = 3, work_dir = tempfile()
  )$traces
  n <- 24
  # COI moves by at most one strain per iteration, so iteration 1 shows the start
  first <- lapply(tr, function(x) as.numeric(x[1, 2:(n + 1)]))
  expect_true(all(first[[1]] <= 2))
  expect_true(all(abs(first[[2]] - 8) <= 1))
  expect_true(all(abs(first[[3]] - 15) <= 1))
})
