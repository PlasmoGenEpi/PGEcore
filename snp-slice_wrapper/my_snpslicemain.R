## execution file for the slice sampler
## this file can
## 1) read input data;
## 2) read command line tuning parameters
## 3) run the mcmc algorithm
## 4) save output to txt;
## last update: 5/15/2023
## ???: nianqiao@purdue.edu

## set the correct working directory before using the code
# setwd()
#
## ## Import necessary libraries
library(tools)
rm(list = ls())

## Parse command-line arguments
args <- commandArgs(trailingOnly = TRUE)

## Set default values for user-determined parameters
input <- "inputdata/example_cat.txt"
input2 <- NA
output_dir <- "output"
model <- 0 ## cat
## prior parameters
rho <- 0.5
alpha <- 2.6
## mcmc parameters
nmcmc <- 10000
gap <- 1
## terminate the chain if lpost does not change in > gap steps
burnin <- 0
## this is the minimum number of MCMC iterations, default is nmcmc/2; see line 47.
rep <- 1
## decide which information to store
store_mcmc <- FALSE
store_sink <- FALSE
store_map <- TRUE

## Parse command-line arguments
for (arg in args) {
  parts <- strsplit(arg, "=")[[1]]
  if (length(parts) == 2) {
    name <- parts[1]
    value <- parts[2]
    if (name == "model") model <- as.integer(value)
    if (name == "nmcmc") nmcmc <- as.integer(value)
    if (name == "alpha") alpha <- as.numeric(value)
    if (name == "input") input <- value
    if (name == "output_dir") output_dir <- value
    if (name == "input2") input2 <- value
    if (name == "gap") gap <- as.integer(value)
  }
}
input_file_name <- tools::file_path_sans_ext(basename(input))
if (is.na(burnin)) burnin <- floor(nmcmc / 2) ## default for burnin is half of nmcmc
## threshold to determine single infections
threshold <- 0.001
## parameters for the observation model
e1 <- 0.05
e2 <- 0.05

set.seed(rep)
## data processing-----
model <- c("cat", "pois", "bin", "neg", "jpois")[model + 1]
if (model == "cat") {
  dfy <- read.delim(input)
  # remove individual id
  dfy$host <- NULL
  dfy$ind_name <- NULL
  N <- dim(dfy)[1]
  P <- dim(dfy)[2]
  y <- as.matrix(dfy)
  r <- NA
} else {
  ## for pois and bin models
  y <- as.matrix(read.delim(input)[, -1])
  r <- as.matrix(read.delim(input2)[, -1]) + y
  N <- dim(y)[1]
  P <- dim(y)[2]
  rho <- sum(y, na.rm = T) / sum(r, na.rm = T)
}
# load helper functions for slice sampler-----
source(paste("source/snpslice", model, ".R", sep = ""))
source("source/snpslicecommon.R")


cat("running", model, "model on data", input_file_name, "\n")
cat("N =", N, "P=", P, "\n")
cat("initialization begin\n")
state <- snp_init(thres = threshold)

cat("initialization done, at", state$loglik, "\n")
cat("starting with", sum(rowSums(state$A) == 1), "single infections\n")

state <- slice_init(state)
cat("plan to run", nmcmc, "iterations of MCMC with", burnin, "iterations for burn in, and gap =", gap, "\n")

map <- list()
lpostmax <- -Inf
mapiter <- 0
mapktrunc <- state$ktrunc
if (store_sink) {
  filename <- paste("sink_", input_file_name, "_", model, "_nmcmc", nmcmc, "_gap", gap, "_rep", rep, ".txt", sep = "")
  sink(file = file.path(output_dir, filename), split = TRUE, type = "output")
}
for (iter in 1:nmcmc) {
  state <- sliceIter(state)
  cat(iter, sum(colSums(state$A) > 0), state$kstar, dim(state$A)[2], sum(rowSums(state$A) == 1), state$logpost, lpostmax, "\n")
  if (state$ktrunc > mapktrunc) { ## if ktrunc changes, restart map
    map <- state
    mapiter <- iter
    mapktrunc <- state$ktrunc
    lpostmax <- state$logpost
  } else if (state$logpost > lpostmax) { ## if still the same ktrunc, but lpost increases,
    map <- state
    mapiter <- iter
    lpostmax <- state$logpost
  }
  if (iter > burnin & mapiter < iter - gap) break
}
if (store_sink) {
  sink()
}

## save the final sample (MAP estimator)
A <- map$A[, which(colSums(map$A) > 0)]
filename <- paste(model, "_A_nmcmc", nmcmc, "_gap", gap, "_rep", rep, ".txt", sep = "")
write.table(A, file.path(output_dir, filename),
  append = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE
)
D <- map$D[which(colSums(map$A) > 0), ]
filename <- paste(model, "_D_nmcmc", nmcmc, "_gap", gap, "_rep", rep, ".txt", sep = "")
write.table(D, file.path(output_dir, filename),
  append = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE
)
## save MCMC state to RData file (for warm starts)
if (store_mcmc) {
  filename <- paste("mcmcRData/mcmc_", "_", model, "_nmcmc", nmcmc, "_gap", gap, "_rep", rep, ".RData", sep = "")
  save(state, rho, alpha, e1, e2, model, r, y, rep, file = filename)
}
