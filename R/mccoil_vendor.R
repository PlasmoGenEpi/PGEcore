# Vendored THEREALMcCOIL R MCMC drivers (categorical + proportional).
# Source: https://github.com/EPPIcenter/THEREALMcCOIL
# C implementations live in src/ and are registered for .C() at install time.

#' THEREALMcCOIL categorical MCMC
#'
#' Ported from `McCOIL_categorical.R`. Calls compiled `McCOIL_categorical`
#' in the PGEcore shared library. Writes MCMC traces and a summary TSV under
#' `path`. Named `run_mccoil_categorical` so it does not collide with the
#' registered C entry point `McCOIL_categorical`.
#'
#' @param data Numeric matrix of heterozygous/homozygous scores (samples x
#'   sites); missing coded as `-1`.
#' @param maxCOI Upper bound for COI.
#' @param threshold_ind Minimum non-missing sites for a sample.
#' @param threshold_site Minimum non-missing samples for a site.
#' @param totalrun Total MCMC iterations.
#' @param burnin Burn-in iterations discarded when summarising.
#' @param M0 Initial COI.
#' @param e1 Probability of calling homozygous loci heterozygous.
#' @param e2 Probability of calling heterozygous loci homozygous.
#' @param err_method `1`/`2` treat error rates as constants; `3` estimates them.
#' @param path Directory for MCMC output files.
#' @param output Base filename for MCMC traces (`<output>_summary.txt` is also
#'   written).
#'
#' @return `NULL`. Called for the files it writes.
#' @keywords internal
run_mccoil_categorical <- function(data,
                               maxCOI = 25,
                               threshold_ind = 20,
                               threshold_site = 20,
                               totalrun = 10000,
                               burnin = 1000,
                               M0 = 15,
                               e1 = 0.05,
                               e2 = 0.05,
                               err_method = 1,
                               path = getwd(),
                               output = "output.txt") {
  In_ind <- rep(NA, nrow(data))
  In_site <- rep(NA, ncol(data))
  for (i in (1:nrow(data))) {
    if ((length(data[i, ]) - sum(data[i, ] == -1)) >= threshold_ind) {
      In_ind[i] <- TRUE
    } else {
      In_ind[i] <- FALSE
    }
  }
  for (i in (1:ncol(data))) {
    if ((length(data[, i]) - sum(data[, i] == -1)) >= threshold_site) {
      In_site[i] <- TRUE
    } else {
      In_site[i] <- FALSE
    }
  }

  ## remove sites and individuals with too much missing data
  simpleS <- data[In_ind, ]
  simpleS <- simpleS[, In_site]

  ## remove sites with P=0 or 1
  P <- rep(NA, ncol(simpleS))
  for (j in (1:ncol(simpleS))) {
    temp <- simpleS[, j]
    P[j] <- (sum(temp == 1) + 0.5 * sum(temp != 0 &
      temp != 1 & temp != -1)) / sum(temp != -1)
  }
  In <- (P != Inf & P != "NaN" & P != 0 & P != 1)
  simpleS2 <- simpleS[, In]

  select_pos <- colnames(data)[In_site][In]
  select_ind <- rownames(data)[In_ind]

  n <- nrow(simpleS2)
  k <- ncol(simpleS2)
  simpleS2_vec <- as.vector(t(simpleS2))
  P0 <- P[In]
  M0 <- rep(M0, n)

  if ((n > 10 & k > 10)) {
    K <- .C(
      "McCOIL_categorical",
      as.integer(maxCOI),
      as.integer(totalrun),
      as.integer(n),
      as.integer(k),
      as.double(simpleS2_vec),
      as.integer(M0),
      as.double(P0),
      as.double(e1),
      as.double(e2),
      as.character(output),
      as.character(path),
      as.integer(err_method),
      PACKAGE = "PGEcore"
    )
  } else {
    stop(paste("Sample size is too small (n=", n, ", k=", k, ").", sep = ""))
  }

  ## summarize results
  outputMCMC1 <- utils::read.table(paste(path, "/", output, sep = ""), header = FALSE)
  meanM <- as.numeric(round(apply(outputMCMC1[(burnin + 1):totalrun, (1:n) + 1], 2, mean)))
  meanP <- as.numeric(apply(outputMCMC1[(burnin + 1):totalrun, ((1:k) + n + 1)], 2, mean))
  medianM <- as.numeric(apply(outputMCMC1[(burnin + 1):totalrun, (1:n) + 1], 2, stats::median))
  medianP <- as.numeric(apply(outputMCMC1[(burnin + 1):totalrun, ((1:k) + n + 1)], 2, stats::median))
  M975 <- as.numeric(apply(outputMCMC1[(burnin + 1):totalrun, (1:n) + 1], 2, function(x) {
    stats::quantile(x, probs = 0.975)
  }))
  P975 <- as.numeric(apply(outputMCMC1[(burnin + 1):totalrun, ((1:k) + n + 1)], 2, function(x) {
    stats::quantile(x, probs = 0.975)
  }))
  M025 <- as.numeric(apply(outputMCMC1[(burnin + 1):totalrun, (1:n) + 1], 2, function(x) {
    stats::quantile(x, probs = 0.025)
  }))
  P025 <- as.numeric(apply(outputMCMC1[(burnin + 1):totalrun, ((1:k) + n + 1)], 2, function(x) {
    stats::quantile(x, probs = 0.025)
  }))
  sdM <- as.numeric(apply(outputMCMC1[(burnin + 1):totalrun, (1:n) + 1], 2, stats::sd))
  sdP <- as.numeric(apply(outputMCMC1[(burnin + 1):totalrun, ((1:k) + n + 1)], 2, stats::sd))

  if (err_method == 3) {
    mean_e1 <- as.numeric(mean(outputMCMC1[(burnin + 1):totalrun, (k + n + 2)]))
    median_e1 <- as.numeric(stats::median(outputMCMC1[(burnin + 1):totalrun, (k + n + 2)]))
    e1_975 <- as.numeric(stats::quantile(outputMCMC1[(burnin + 1):totalrun, (k + n + 2)], probs = 0.975))
    e1_025 <- as.numeric(stats::quantile(outputMCMC1[(burnin + 1):totalrun, (k + n + 2)], probs = 0.025))
    sd_e1 <- as.numeric(stats::sd(outputMCMC1[(burnin + 1):totalrun, (k + n + 2)]))

    mean_e2 <- as.numeric(mean(outputMCMC1[(burnin + 1):totalrun, (k + n + 3)]))
    median_e2 <- as.numeric(stats::median(outputMCMC1[(burnin + 1):totalrun, (k + n + 3)]))
    e2_975 <- as.numeric(stats::quantile(outputMCMC1[(burnin + 1):totalrun, (k + n + 3)], probs = 0.975))
    e2_025 <- as.numeric(stats::quantile(outputMCMC1[(burnin + 1):totalrun, (k + n + 3)], probs = 0.025))
    sd_e2 <- as.numeric(stats::sd(outputMCMC1[(burnin + 1):totalrun, (k + n + 3)]))
  }
  if ((err_method == 1) | (err_method == 2)) {
    output_sum <- data.frame(cbind(
      rep(output, (n + k)),
      c(rep("C", n), rep("P", k)),
      c(select_ind, select_pos),
      c(meanM, meanP),
      c(medianM, medianP),
      round(c(sdM, sdP), digits = 5),
      c(M025, P025),
      c(M975, P975)
    ))
  } else {
    output_sum <- data.frame(cbind(
      rep(output, (n + k + 2)),
      c(rep("C", n), rep("P", k), "e1", "e2"),
      c(select_ind, select_pos, "e1", "e2"),
      c(meanM, meanP, mean_e1, mean_e2),
      c(medianM, medianP, median_e1, median_e2),
      round(c(sdM, sdP, sd_e1, sd_e2), digits = 5),
      c(M025, P025, e1_025, e2_025),
      c(M975, P975, e1_975, e2_975)
    ))
  }
  colnames(output_sum) <- c(
    "file",
    "CorP",
    "name",
    "mean",
    "median",
    "sd",
    "quantile0.025",
    "quantile0.975"
  )
  utils::write.table(
    output_sum,
    paste(path, "/", output, "_summary.txt", sep = ""),
    sep = "\t",
    col.names = TRUE,
    row.names = FALSE,
    quote = FALSE
  )
  invisible(K)
}

#' THEREALMcCOIL proportional MCMC
#'
#' Ported from `McCOIL_proportional.R`. Reads the fitted beta grid from package
#' `extdata` and calls compiled `McCOIL_prop`. Named `run_mccoil_proportional`
#' so it does not collide with registered C symbols.
#'
#' @param dataA1 Allele-1 read counts (samples x sites).
#' @param dataA2 Allele-2 read counts (samples x sites).
#' @param maxCOI Upper bound for COI.
#' @param totalrun Total MCMC iterations.
#' @param burnin Burn-in iterations discarded when summarising.
#' @param M0 Initial COI.
#' @param epsilon Sequencing error parameter for the proportional model.
#' @param err_method `1`/`2` treat epsilon as constant; `3` estimates it.
#' @param path Directory for MCMC output files.
#' @param output Base filename for MCMC traces.
#'
#' @return `NULL`. Called for the files it writes.
#' @keywords internal
run_mccoil_proportional <- function(dataA1,
                                dataA2,
                                maxCOI = 25,
                                totalrun = 10000,
                                burnin = 1000,
                                M0 = 15,
                                epsilon = 0.02,
                                err_method = 1,
                                path = getwd(),
                                output = "output.txt") {
  grid_path <- system.file(
    "extdata",
    "fitted_beta_grid_25.txt",
    package = "PGEcore",
    mustWork = TRUE
  )
  grid <- utils::read.table(grid_path, header = TRUE)
  n <- nrow(dataA1)
  k <- ncol(dataA1)
  M0 <- rep(M0, n)
  P0 <- rep(0.5, k)
  A1 <- as.vector(t(dataA1))
  A2 <- as.vector(t(dataA2))

  if ((n > 10 & k > 10)) {
    Kc <- .C(
      "McCOIL_prop",
      as.integer(maxCOI),
      as.integer(totalrun),
      as.integer(n),
      as.integer(k),
      as.double(A1),
      as.double(A2),
      as.integer(M0),
      as.double(P0),
      as.double(grid$A),
      as.double(grid$B),
      as.double(epsilon),
      as.character(output),
      as.character(path),
      as.integer(err_method),
      PACKAGE = "PGEcore"
    )
  } else {
    stop(paste("Sample size is too small (n=", n, ", k=", k, ").", sep = ""))
  }

  ## summarize results
  outputMCMC2 <- utils::read.table(paste(path, "/", output, sep = ""), header = FALSE)
  meanM <- as.numeric(round(apply(outputMCMC2[(burnin + 1):totalrun, (1:n) + 1], 2, mean)))
  meanP <- as.numeric(apply(outputMCMC2[(burnin + 1):totalrun, ((1:k) + n + 1)], 2, mean))
  medianM <- as.numeric(apply(outputMCMC2[(burnin + 1):totalrun, (1:n) + 1], 2, stats::median))
  medianP <- as.numeric(apply(outputMCMC2[(burnin + 1):totalrun, ((1:k) + n + 1)], 2, stats::median))
  M975 <- as.numeric(apply(outputMCMC2[(burnin + 1):totalrun, (1:n) + 1], 2, function(x) {
    stats::quantile(x, probs = 0.975)
  }))
  P975 <- as.numeric(apply(outputMCMC2[(burnin + 1):totalrun, ((1:k) + n + 1)], 2, function(x) {
    stats::quantile(x, probs = 0.975)
  }))
  M025 <- as.numeric(apply(outputMCMC2[(burnin + 1):totalrun, (1:n) + 1], 2, function(x) {
    stats::quantile(x, probs = 0.025)
  }))
  P025 <- as.numeric(apply(outputMCMC2[(burnin + 1):totalrun, ((1:k) + n + 1)], 2, function(x) {
    stats::quantile(x, probs = 0.025)
  }))
  sdM <- as.numeric(apply(outputMCMC2[(burnin + 1):totalrun, (1:n) + 1], 2, stats::sd))
  sdP <- as.numeric(apply(outputMCMC2[(burnin + 1):totalrun, ((1:k) + n + 1)], 2, stats::sd))

  if (err_method == 3) {
    mean_e3 <- as.numeric(mean(outputMCMC2[(burnin + 1):totalrun, (k + n + 2)]))
    median_e3 <- as.numeric(stats::median(outputMCMC2[(burnin + 1):totalrun, (k + n + 2)]))
    e3_975 <- as.numeric(stats::quantile(outputMCMC2[(burnin + 1):totalrun, (k + n + 2)], probs = 0.975))
    e3_025 <- as.numeric(stats::quantile(outputMCMC2[(burnin + 1):totalrun, (k + n + 2)], probs = 0.025))
    sd_e3 <- as.numeric(stats::sd(outputMCMC2[(burnin + 1):totalrun, (k + n + 2)]))
  }
  if ((err_method == 1) | (err_method == 2)) {
    output_sum <- data.frame(cbind(
      rep(output, (n + k)),
      c(rep("C", n), rep("P", k)),
      c(rownames(dataA1), colnames(dataA1)),
      c(meanM, meanP),
      c(medianM, medianP),
      round(c(sdM, sdP), digits = 5),
      c(M025, P025),
      c(M975, P975)
    ))
  } else {
    output_sum <- data.frame(cbind(
      rep(output, (n + k + 1)),
      c(rep("C", n), rep("P", k), "epsilon"),
      c(rownames(dataA1), colnames(dataA1), "epsilon"),
      c(meanM, meanP, mean_e3),
      c(medianM, medianP, median_e3),
      round(c(sdM, sdP, sd_e3), digits = 5),
      c(M025, P025, e3_025),
      c(M975, P975, e3_975)
    ))
  }
  colnames(output_sum) <- c(
    "file", "CorP", "name", "mean", "median", "sd",
    "quantile0.025", "quantile0.975"
  )
  utils::write.table(
    output_sum,
    paste(path, "/", output, "_summary.txt", sep = ""),
    sep = "\t",
    col.names = TRUE,
    row.names = FALSE,
    quote = FALSE
  )
  invisible(Kc)
}
