library(malaria.em)
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
    "--subset_alleles",
    help = "A logical to indicate whether the matrix should be subset by user-supplied allele_names",
    type = "logical",
    default = FALSE
  ),
  make_option(
    "--allele_names",
    help = "Comma-separated allele names to subset (e.g. 'msp1,ama1,csp'). Use 'NULL' to skip subsetting.",
    type = "character",
    default = NULL
  ),
  make_option(
    "--coi",
    help = "Single COI value (e.g. '2') or comma-separated values (e.g. '1,2,3') for a COI distribution. This should be the largest COI observed in the data. If length is greater than 1, a zero truncated Poisson distribution on COI is assumed",
    type = "character",
    default = "1, 2, 3"
  ),
  make_option(
    "--write_results",
    help = str_c("A logical to indicate whether the malaria.em results should be written to file (.tsv). If TRUE ",
                 # Population-level multilocus genotype frequency estimate + standard error 
                 "- 'gt_freq_summary.tsv', with columns: ",
                 "gt_id, target_id, seq, freq, freq_se",
                 # Sample-level phased genotypes + posterior probability
                 "- 'gt_phase_summary.tsv', with columns: ",
                 "specimen_id, target_id, seq, gt_id, posterior_est, phase_id"),
    type = "logical",
    default = TRUE
  )
)

arg <- parse_args(OptionParser(option_list = opts))

# Convert allele_names from string to character vector or NULL
allele_names <- if (arg$allele_names %in% c("NULL", "", NULL)) {
  NULL
} else {
  strsplit(arg$allele_names, ",")[[1]]
}

# Convert coi from string to numeric vector
coi <- as.numeric(strsplit(arg$coi, ",")[[1]])

required_args <- c("allele_table", "coi")
missing_args <- setdiff(required_args, names(arg))
if (length(missing_args) > 0) {
  stop(stringr::str_c("Missing required arguments: ", paste(missing_args, collapse = ", ")))
}


# {malaria.em} input function ---------------------------------------------

#' Load and format allele table for malaria.em
#' 
#' This function loads the allele table in .tsv format and reshapes it into a matrix suitable as input for the `malaria.em` package.
#' The input file should contain columns named `specimen_id`, `target_id`, and `seq`. The function pivots the data to wide format
#' and optionally filters to retain only a specified subset of alleles (e.g., the most diverse loci).
#'
#' @param file_path Character. Path to the input .tsv file containing allele calls with columns `specimen_id`, `target_id`, and `seq`. 
#' @param subset_alleles Logical. Whether to subset the allele table to specific loci/alleles (default is `FALSE`).
#' @param allele_names Character vector. Optional vector of allele names (column names) to retain if `subset_alleles = TRUE`.
#'
#' @return A matrix where rows represent specimens/samples and columns represent alleles, with each cell containing the unique sequence.
#' This matrix is formatted for use with `malaria.em::malaria_em()` and wrapper function `run_malariaem()`.
#' 
#' @importFrom readr read_tsv cols col_character
#' @importFrom dplyr select group_by any_of
#' @importFrom tidyr pivot_wider
#' @importFrom tibble column_to_rownames
#' @importFrom checkmate assert_file_exists assert_character
#' 
#' @export
#'
#' @examples
#' \dontrun{
#' # Load all alleles
#' matrix <- load_allele_table("data/example2_allele_table.tsv")
#'
#' # Load only a subset of alleles (e.g., most diverse loci)
#' top_loci <- c("Pf3D7_12_v3-0659902-0660096", "Pf3D7_13_v3-1419395-1419601")
#' matrix_subset <- load_allele_table("data/example2_allele_table.tsv", subset_alleles = TRUE, allele_names = top_loci)
#' }
load_allele_table <- function(file_path, subset_alleles = FALSE, allele_names){
  
  # Check inputs
  checkmate::assert_file_exists(file_path, access = "r")
  if (subset_alleles) {
    checkmate::assert_character(allele_names, any.missing = FALSE)
  }

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
  
  mhap_wide <- mhap |>
    select(specimen_id, target_id, seq) |>
    group_by(specimen_id, target_id) |>
    pivot_wider(names_from = target_id,
                values_from = seq,
                values_fn = ~ paste(unique(.), collapse = " ")) |>
    column_to_rownames(var = "specimen_id")
  
  if (subset_alleles) {
    print("Subsetting alleles...")
    non_existing <- setdiff(allele_names, colnames(mhap_wide))
    if (length(non_existing) > 0) {
      stop("The following allele_names were not found in the data: ", paste(non_existing, collapse = ", "))
    }
    mhap_final <- mhap_wide |> dplyr::select(dplyr::any_of(allele_names))
  } else {
    mhap_final <- mhap_wide
  }
  
  # if(subset_alleles == TRUE){
  # 
  #   print("Subsetting alleles...")
  # 
  #   mhap_final <- mhap_wide |>
  #       select(any_of(allele_names))
  # }
  # 
  # else{
  #   mhap_final <- mhap_wide
  # }
  
  print("Converting to matrix...")
    
  mat <- mhap_final |>
    as.matrix()
  
  return(mat)
}


# malaria.em run function -------------------------------------------------

#' Run malaria.em and summarize results
#' 
#' Executes the `malaria.em()` function on an allele matrix and a vector of complexity of infection (COI) values.
#' It returns the default `malaria.em()` output object as well as summarized outputs for population-level multi-locus genotype frequencies and sample-level multi-locus genotype phase estimates. It also optionally writes the summary outputs to a .tsv file.
#' *NOTE: the package can become very slow if there are many possible haplotype combinations. It is recommended to start by running a subset of alleles to gauge computational requirements.*
#'
#' @param matrix A matrix of allele sequences (specimens as rows, loci as columns), formatted as required by `malaria.em`. This can be generated by using the wrapper function `load_allele_table()`
#' @param coi An integer or numeric vector of possible COI (complexity of infection) range. If a vector is provided, the estimation will assume a zero-truncated Poisson distribution on COI. Otherwise, the estimation will assume COI is fixed.
#' @param write_results Logical. Whether to write the genotype frequency and phase summary results to a .tsv file (default is `TRUE`).
#' 
#' #' @details
#' If \code{write_results = TRUE}, the function will write two tab-separated files to the current working directory:
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
#' 
#' @return A list containing three objects:
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
#' }
#' @importFrom dplyr select distinct left_join rename mutate
#' @importFrom malaria.em malaria.em
#' @importFrom tibble rowid_to_column rownames_to_column
#' @importFrom tidyr pivot_longer separate_rows
#' @importFrom readr write_tsv
#' @importFrom checkmate assert_matrix assert_numeric assert_logical
#' 
#' @export
#'
#' @examples
#' \dontrun{
#' # Example matrix and COI
#' matrix <- load_allele_table("data/example2_allele_table.tsv")
#' coi <- c(1:4) # assume COI ranges from 1 to 4
#' result <- run_malariaem(matrix, coi, write_results = FALSE)
#' }
#' 

run_malariaem <- function(matrix, coi, write_results = TRUE){
  
  # Checks inputs
  checkmate::assert_matrix(matrix)
  checkmate::assert_numeric(coi, any.missing = FALSE)

  sample_name <- matrix |> 
    as.data.frame() |> 
    rownames_to_column("specimen_id") |> 
    distinct(specimen_id) |> 
    rowid_to_column("ids") 
  
  allele_names <- colnames(matrix)
  
  print("Running malaria.em...")
  output <- malaria.em::malaria.em(matrix, sizes = coi, locus.label = allele_names)

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
  
  
  print("Summarizing sample-level phased multi-locus genotypes + posterior probability estimates...")
  # Phased genotypes per sample + posterior probability
  haplo_pred <- as.data.frame(output$pred.haplo.set) |> rename(haplo.set = `output$pred.haplo.set`) |> rowid_to_column("ids")
  haplo_pred_prob <- as.data.frame(output$haplo.sets) 
  
  haplo_wide <- as.data.frame(output$haplotype) |> rowid_to_column("gt_id") 
  
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
  
  if(write_results){
  print("Writing results to file...")
  readr::write_tsv(gt_freq_summary, "gt_freq_summary.tsv")
  readr::write_tsv(gt_phase_summary, "gt_phase_summary.tsv")
  }
  else{
    NULL
  }
  
  return(list(output, gt_freq_summary, gt_phase_summary))
  
}

# MAIN --------------------------------------------------------------------

matrix <- load_allele_table(file_path = arg$allele_table, 
                            subset_alleles = arg$subset_alleles,
                            allele_names = allele_names)
results <- run_malariaem(matrix = matrix,
                         coi = coi, 
                         write_results = arg$write_results)


# TEST --------------------------------------------------------------------
if(interactive()){

alleles <- c("Pf3D7_12_v3-0659902-0660096", "Pf3D7_13_v3-1419395-1419601")
coi = c(1:4)

matrix <- load_allele_table(file_path = "data/example2_allele_table.tsv",
                            subset_alleles = TRUE,
                            allele_names = alleles)
out <- run_malariaem(matrix, coi, write_results = TRUE)

}