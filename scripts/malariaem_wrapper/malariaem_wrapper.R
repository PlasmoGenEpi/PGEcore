library(malaria.em)
library(readr)
library(dplyr)
library(tidyr)
library(tibble)
library(stringr)
library(optparse)
library(validate)
library(checkmate)

# Parse arguments ---------------------------------------------------------

options(dplyr.summarise.inform = FALSE)

opts <- list(
  make_option(
    "--allele_table",
    help = str_c("TSV containing alleles present in each specimen, with columns: ",
                 "specimen_id, target_id, seq"),
    type = "character",
    default = NULL,
    callback = function(opt, flag_string, value, parser, ...){
      if (!file.exists(value)) {
        stop(stringr::str_c(value, " does not exist"))
      }
      value
    }
  ),
  make_option(
    "--subset_targets",
    help = "A logical to indicate whether the matrix should be subset by user-supplied target_id",
    type = "logical",
    default = FALSE
  ),
  make_option(
    "--target_groups",
    help = str_c("TSV containing the targets that should be analyzed for specific groups (eg pfdhr/pfdhps), with columns: ",
                 "group_id, target_id"),
    type = "character",
    default = NULL,
    callback = function(opt, flag_string, value, parser, ...){
      if (!file.exists(value)) {
        stop(stringr::str_c(value, " does not exist"))
      }
      value
    }
  ),
  make_option(
    "--max_size",
    help = str_c("Single population-COI value (e.g. '2')",
                 "Cutoff to not run Malaria EM due to time constraints",
                 "Default is max size of 8"),
    type = "character",
    default = "8"
  ),
  make_option(
    "--test_size",
    help = str_c("Single COI value (e.g. '2'), or 'min' to default to the minimum allowed set size",
                 "Default is 'min'"),
    type = "character",
    default = "min"
  ),
  make_option(
    "--freq_output",
    help = "Path to frequency output file",
    type = "character",
    default = "gt_freq_summary_all.tsv",
  ),
  make_option(
    "--phase_output",
    help = "Path to phasing output file",
    type = "character",
    default = "gt_phase_summary_all.tsv",
  ),
  make_option(
    "--seed",
    help = "seed",
    type = "character",
    default = "1",
  )
)

arg <- parse_args(OptionParser(option_list = opts))
set.seed(as.numeric(arg$seed))

required_args <- c("allele_table")
missing_args <- setdiff(required_args, names(arg))
if (length(missing_args) > 0) {
  stop(stringr::str_c("Missing required arguments: ", paste(missing_args, collapse = ", ")))
}

# {malaria.em} input functions ---------------------------------------------

#' Load and format allele table for malaria.em
#' 
#' This function loads the allele table in .tsv format and reshapes it into a matrix suitable as input for the `malaria.em` package.
#' The input file should contain columns named `specimen_id`, `target_id`, and `seq`. 
#'
#' @param file_path Character. Path to the input .tsv file containing allele calls with columns `specimen_id`, `target_id`, and `seq`. 
#'
#' @return A matrix where rows represent specimens/samples and columns represent alleles, with each cell containing the unique sequence.
#' This matrix is formatted for use with `malaria.em::malaria_em()` and wrapper function `run_malariaem()`.
#' 
#' @importFrom readr read_tsv cols col_character
#' @importFrom dplyr select group_by
#' @importFrom tidyr pivot_wider
#' @importFrom tibble column_to_rownames
#' @importFrom checkmate assert_file_exists assert_character
#' 
#' @export
#'
#' @examples
#' \dontrun{
#' # Load allele table
#' matrix <- load_allele_table("data/example2_allele_table.tsv")
#' 
#' }
load_allele_table <- function(file_path){
  
  # Check inputs
  checkmate::assert_file_exists(file_path, access = "r")

  print("Reading input data...")
  
  mhap <- readr::read_tsv(
    file_path,
    col_types = readr::cols(
      specimen_id = readr::col_character(),
      target_id = readr::col_character(),
      seq = readr::col_character()
    ),
    col_select = c("specimen_id", "target_id", "seq")
  )
  
  # Validate input data
  rules <- validate::validator(
    is.character(specimen_id), 
    is.character(target_id), 
    is.character(seq), 
    ! is.na(specimen_id), 
    ! is.na(target_id), 
    ! is.na(seq)
  )
  fails <- validate::confront(mhap, rules, raise = "all") %>%
    validate::summary() %>%
    dplyr::filter(fails > 0)
  if (nrow(fails) > 0) {
    stop(
      "Input input_data failed one or more validation checks: ", 
      str_c(fails$expression, collapse = "\n"), 
      call. = FALSE
    )
  }
  
  print("Converting to matrix...")
  
  matrix <- mhap |>
    select(specimen_id, target_id, seq) |>
    group_by(specimen_id, target_id) |>
    pivot_wider(names_from = target_id,
                values_from = seq,
                values_fn = ~ paste(unique(.), collapse = " ")) |>
    column_to_rownames(var = "specimen_id") |>
    as.matrix()
  
  return(matrix)
}

#' Load target groups data
#' 
#' This function loads the target groups data in .tsv format. The input file should contain columns named `group_id` and `target_id`.
#'
#' @param file_path Character. Path to the input .tsv file with columns `group_id` and `target_id`.
#'
#' @return A dataframe containing groups of targets to be analyzed separately by `malaria.em`.
#' This object can be optionally supplied to the wrapper function `run_malariaem()` to run `malaria.em` on subsets of targets. 
#' 
#' @importFrom readr read_tsv cols col_character
#' @importFrom checkmate assert_file_exists assert_character
#' @importFrom validate validator confront summary
#' @importFrom stringr str_c
#' 
#' @export
#' 
#' @examples
#' \dontrun{
#' # Load allele table
#' target_groups <- load_target_groups("data/example_target_groups.tsv")
#' }
load_target_groups <- function(file_path) {
  
  # Check inputs
  checkmate::assert_file_exists(file_path, access = "r")
  
  print("Reading target groups...")
  
  target_groups <- read_tsv(file_path, 
                            col_types = cols(group_id = col_character(),
                                             target_id = col_character()
                            ),
                            show_col_types = FALSE
  )
  
  # Validate
  rules <- validate::validator(
    is.character(group_id), 
    is.character(target_id), 
    ! is.na(group_id), 
    ! is.na(target_id) 
  )
  fails <- validate::confront(target_groups, rules, raise = "all") %>%
    validate::summary() %>%
    dplyr::filter(fails > 0)
  if (nrow(fails) > 0) {
    stop(
      "Input input_data failed one or more validation checks: ", 
      str_c(fails$expression, collapse = "\n"), 
      call. = FALSE
    )
  }
  
  return(target_groups)
}

# malaria.em run function -------------------------------------------------

#' Run malaria.em and summarize results
#' 
#' Executes the `malaria.em()` function on an allele matrix and a vector (or single value) of complexity of infection (COI).
#' It returns the default `malaria.em()` output object as well as summarized outputs for population-level multi-locus genotype frequencies and sample-level multi-locus genotype phase estimates. 
#' It also writes two TSV summary outputs to `output_dir` (or the current working directory if `output_dir` is not specified).
#' *NOTE: the package can become very slow if there are many possible haplotype combinations. It is recommended to start by running a subset of alleles to gauge computational requirements. This can be specified through by supplying `subset_targets = TRUE` and `target_groups`*
#'
#' @param matrix A matrix of allele sequences (specimens as rows, loci as columns), formatted as required by `malaria.em`. This can be generated by using the wrapper function `load_allele_table()`
#' @param test_size Value of the maximum size to test or "min" to use the minimum allowed level. Overridden by the minimum allowed.
#' @param max_size Value of the maximum size to throw an error in case of an override.
#' @param subset_targets Logical. If `TRUE`, run `malaria.em` separately for each group defined in `target_groups`, otherwise run once on all columns in the `matrix`. Default is `FALSE`.
#' @param target_groups A data frame with columns `group_id` and `target_id`. When `subset_targets = TRUE`, `target_groups` must be supplied and each unique `group_id` defines a set of targets to include in a separate run. All `target_id` values must be column names in `matrix`.
#' @param output_dir Character path to the directory where TSV outputs are written (note this directory must already exist). If `NULL`, uses `getwd()` and outputs to the current working directory. 
#' 
#' #' @details
#' For each run (either single full run, or one per `group_id` when `subset_targets = TRUE`), the function will write two tab-separated files to the current working directory or `output_dir`:
#' \itemize{
#'   \item \code{gt_freq_summary.tsv}: Population-level multilocus genotype frequency estimates. Columns:
#'     \itemize{
#'       \item \code{gt_id}: Unique genotype ID (global across all samples).
#'       \item \code{target_id}: Allele name.
#'       \item \code{seq}: Allele sequence at the locus.
#'       \item \code{freq}: Estimated population-level frequency of the multilocus genotype.
#'       \item \code{freq_se}: Standard error of the estimated frequency.
#'     }
#'   \item \code{gt_phase_summary.tsv}: Sample-level phased multilocus genotypes with posterior probabilities. Columns:
#'     \itemize{
#'       \item \code{specimen_id}: Sample identifier.
#'       \item \code{target_id}: Allele name.
#'       \item \code{seq}: Allele sequence at the locus.
#'       \item \code{gt_id}: Multi-locus genotype ID (global across samples). 
#'       \item \code{posterior_est}: Posterior probability of the phased genotype.
#'       \item \code{phase_id}: A sample-specific numeric label (1, 2, ...) assigned to each phased multi-locus genotype within a sample. This is independent of the global \code{gt_id} and is used only for local reference.
#'     }
#' }
#' When subsetting, filenames include a suffix of the `group_id` (non-alphanumeric characters replaced by underscores).
#' 
#' @return If `subset_targets = FALSE` a single list containing three objects:
#' \itemize{
#'   \item \code{output}: The full result object returned by `malaria.em::malaria.em()`, which includes:
#'   \itemize{
#'     \item \code{haplo.prob.tab}: A matrix of unique haplotypes with estimated probabilities and standard errors.
#'     \item \code{haplotype}: A matrix where each row is a unique haplotype (columns = loci).
#'     \item \code{haplo.prob}: A vector of MLEs (maximum likelihood estimates) of haplotype probabilities.
#'     \item \code{haplo.prob.std}: A vector of standard errors corresponding to each estimated haplotype frequency.
#'     \item \code{lambda}: Estimated Poisson parameter (if multiple COI values are provided).
#'     \item \code{NumofInfection}: Estimated number of infections.
#'     \item \code{haplo.sets}: A data frame of all possible haplotype combinations and posterior probabilities per subject.
#'     \item \code{n.haplo.set}: A vector indicating the number of possible haplotype sets per subject.
#'     \item \code{pred.haplo.set}: A list of predicted haplotype combinations (by row index in `haplotype`) per subject.
#'   }
#'   \item \code{gt_freq_summary}: A data frame of population-level multilocus genotype frequencies and standard errors.
#'   \item \code{gt_phase_summary}: A data frame summarizing phased multilocus genotypes per sample with posterior estimates. Outputs only one phased MLG per sample (i.e. not the full per-sample MLG posterior probability distribution)
#'   \item \code{group_id}: `NULL` for the full run (no subsetting).
#' }
#' If \code{subset_targets = TRUE}, a named list where each element corresponds to one
#' \code{group_id} and contains the same fields as above (with \code{group_id} set to that label).
#' 
#' @importFrom dplyr select distinct left_join rename mutate dense_rank join_by
#' @importFrom malaria.em malaria.em
#' @importFrom tibble rowid_to_column rownames_to_column
#' @importFrom tidyr pivot_longer separate_rows
#' @importFrom readr write_tsv
#' @importFrom checkmate assert_matrix assert_numeric assert_logical
#' @importFrom stats setNames
#' 
#' @seealso \code{\link{load_allele_table}}, \code{\link[malaria.em]{malaria.em}}
#' 
#' @export
#'
#' @examples
#' \dontrun{
#' # Full run on all targets; write results to a specific directory
#' mat <- load_allele_table("data/example2_allele_table.tsv")
#' coi_range <- c(1:4) # assume COI ranges from 1 to 4
#' result <- run_malariaem(
#'   matrix = mat,
#'   coi_range = coi_range,
#'   subset_targets = FALSE,
#'   output_dir = "results/malariaem_full"
#' )
#' 
#' # Full run on all targets; write results to current working directory
#' mat <- load_allele_table("data/example2_allele_table.tsv")
#' coi_range <- c(1:4) # assume COI ranges from 1 to 4
#' result <- run_malariaem(
#'   matrix = mat,
#'   coi_range = coi_range,
#'   subset_targets = FALSE
#' )
#'
#' # Subset by groups of targets (two groups) and write one set per group
#' mat <- load_allele_table("data/example2_allele_table.tsv")
#' groups <- load_target_groups("data/example_target_groups.tsv)
#' coi_range <- c(1:4) # assume COI ranges from 1 to 4
#' 
#' result_list <- run_malariaem(
#'   matrix = mat,
#'   coi_range = coi_range,
#'   subset_targets = TRUE,
#'   target_groups = groups,
#'   output_dir = "results/malariaem_groups"
#' )
#' }
#' 

run_malariaem <- function(matrix, test_size, max_size, subset_targets = FALSE, target_groups = NULL, freq_output = NULL, phase_output = NULL){
  
  # Checks inputs
  checkmate::assert_matrix(matrix)
  checkmate::assert_flag(subset_targets)
  if (subset_targets) {
    checkmate::assert_data_frame(target_groups, min.rows = 1, min.cols = 2)
    checkmate::assert_subset(c("group_id", "target_id"), choices = names(target_groups))
  }
  
  # Set output directory for writing results
  freq_outdir <- if(is.null(freq_output)) getwd() else dirname(freq_output)
  phase_outdir <- if(is.null(phase_output)) getwd() else dirname(phase_output)

  # Create output dirs if they don't exist
  if(!dir.exists(freq_outdir)) dir.create(freq_outdir, recursive = TRUE, showWarnings = FALSE)
  if(!dir.exists(phase_outdir)) dir.create(phase_outdir, recursive = TRUE, showWarnings = FALSE)

  checkmate::assert_directory_exists(freq_outdir, access = "w")
  checkmate::assert_directory_exists(phase_outdir, access = "w")
  
  run_malariaem_helper <- function(matrix, coi_range, label = NULL){
    min_targets <- max(apply(matrix, 2,num.unique.allele.locus))
    if(test_size == "min"){
      coi_range <- c(1:min_targets)
    }
    else{
      top <- as.numeric(test_size)
      if(min_targets > top){
        top <- min_targets
      }
      coi_range <- c(1:top)}
    
    if(max(coi_range) > max_size){
      stop(
        "Too many alleles to run EM efficiently", 
        call. = FALSE
      )
    }
    #removes samples that are missing at any loci
    matrix <- matrix[rowSums(is.na(matrix)) == 0,]
    
    # get sample names
    sample_name <- matrix |> 
      as.data.frame() |> 
      rownames_to_column("specimen_id") |> 
      distinct(specimen_id) |> 
      rowid_to_column("ids") 
    
    # get target names
    target_names <- colnames(matrix)
    
    message("Running malaria.em", if (!is.null(label)) paste0(" [", label, "]"), "...")
    output <- malaria.em::malaria.em(matrix, sizes = coi_range, locus.label = target_names)
    
    print("Summarizing population-level multi-locus genotype frequency estimates + SE...")
    # Population-level multilocus genotype frequency estimate + standard error 
    gt_freq_summary <- output$haplo.prob.tab |> 
      as.data.frame() |>
      rowid_to_column("gt_id") |>
      pivot_longer(cols = -c("gt_id", "hap.prob", "hap.prob.std"), 
                   names_to = "hap_id", 
                   values_to = "seq") |>
      select(gt_id, 
             target_id = hap_id, 
             seq, 
             freq = hap.prob, 
             freq_se = hap.prob.std)
    
    if (!is.null(label)) {
      gt_freq_summary <- gt_freq_summary |> mutate(group_id = label)
    }
    
    print("Summarizing sample-level phased multi-locus genotypes + posterior probability estimates...")
    # Phased genotypes per sample + posterior probability
    haplo_pred <- as.data.frame(output$pred.haplo.set) |> rename(haplo.set = `output$pred.haplo.set`) |> rowid_to_column("ids")
    haplo_pred_prob <- as.data.frame(output$haplo.sets) 
    
    gt_phase_summary <- haplo_pred |>
      as.data.frame() |> 
      left_join(haplo_pred_prob, 
                by = join_by(ids, haplo.set)) |>
      left_join(sample_name, join_by(ids)) |>
      select(specimen_id,
             haplo.set,
             posterior_est = post.p) |> 
      separate_rows(haplo.set, sep = " ") |> 
      rename(gt_id = haplo.set) |>
      mutate(gt_id = as.integer(gt_id)) |>
      left_join(gt_freq_summary, join_by(gt_id), relationship = "many-to-many") |>
      select(specimen_id, target_id, seq, gt_id, posterior_est) |> 
      mutate(phase_id = dense_rank(gt_id), .by = specimen_id)
    
    if (!is.null(label)) {
      gt_phase_summary <- gt_phase_summary |> mutate(group_id = label)
    }

    list(
      output = output,
      gt_freq_summary = gt_freq_summary,
      gt_phase_summary = gt_phase_summary,
      group_id = label
    )
  }
  
  # Run all with no subsetting
  if(!subset_targets){
    res <- run_malariaem_helper(matrix, coi_range, label = NULL)
    readr::write_tsv(res$gt_freq_summary, freq_output)
    readr::write_tsv(res$gt_phase_summary, phase_output)
    return(res)
  }
  
  # Run with subsetting by group id
  group_ids <- unique(target_groups$group_id) |> as.character()
  
  results_list <- setNames(vector("list", length(group_ids)), group_ids)
  
  for (gid in group_ids) {
    # get targets for this gid
    targets <- unique(target_groups$target_id[target_groups$group_id == gid])
    
    if (length(targets) == 0L) {
      warning("No targets found for group_id '", gid, "'. Skipping.")
      next
    }
    
    # check all exist in matrix
    non_existing <- setdiff(targets, colnames(matrix))
    if (length(non_existing) > 0) {
      stop("Group '", gid, "' requests targets not present in matrix: ",
           paste(non_existing, collapse = ", "))
    }
    
    matrix_subset <- matrix[, targets, drop = FALSE]
    if (ncol(matrix_subset) == 0) {
      warning("Matrix for group_id '", gid, "' is empty; skipping.")
      next
    }
    
    results_list[[gid]] <- run_malariaem_helper(matrix_subset, coi_range, label = gid)
  }
  # Combine results from all runs
  all_gt_freq_summary <- bind_rows(lapply(results_list, `[[`, "gt_freq_summary"))
  all_gt_phase_summary <- bind_rows(lapply(results_list, `[[`, "gt_phase_summary"))
  
  # Write single combined files
  readr::write_tsv(all_gt_freq_summary, freq_output)
  readr::write_tsv(all_gt_phase_summary, phase_output)
  
  return(results_list)
  
}

# MAIN --------------------------------------------------------------------

matrix <- load_allele_table(file_path = arg$allele_table)
if (arg$subset_targets){
  target_groups <- load_target_groups(file_path = arg$target_groups)

  results <- run_malariaem(matrix = matrix,
                          test_size = arg$test_size,
                          max_size = arg$max_size,
                          subset_targets = arg$subset_targets,
                          target_groups,
                          freq_output=arg$freq_output,
                          phase_output=arg$phase_output)
} else {
  results <- run_malariaem(matrix = matrix,
                           test_size = arg$test_size,
                           max_size = arg$max_size,
                          subset_targets = arg$subset_targets,
                          freq_output=arg$freq_output,
                          phase_output=arg$phase_output)
}

# TEST --------------------------------------------------------------------
if(interactive()){

matrix <- load_allele_table(file_path = "data/example2_allele_table.tsv")

target_groups <- load_target_groups(file_path = "data/example_target_groups.tsv")

coi_range <- c(1:4)

out <- run_malariaem(matrix, coi_range, subset_targets = TRUE, target_groups)

}
