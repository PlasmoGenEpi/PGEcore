#!/usr/bin/env Rscript


packagesToLoad = c("tibble", "dplyr", "tidyr", "stringr", "readr", "optparse")

loaded = lapply(packagesToLoad, library, warn.conflicts = F, character.only = TRUE)

options(readr.show_col_types = FALSE)
options(dplyr.summarise.inform = FALSE)


#' Check for required arguments, and report which are missing 
#'
#' @param arg the parsed arguments from arg parse
#' @param required_args the required arguments
#'
#' @return returns void if all required arguments
checkOptparseRequiredArgsThrow <- function(arg, required_args){
  missing <- setdiff(required_args, names(arg))
  if(length(missing) > 0){
    missing = paste0("--", missing)
    stop(paste0("mssing the following arguments: ", paste0(missing, collapse = ", ")))
  }
}


#' Check sharing and unqiue values between two vectors
#'
#' @param vectorA the first vector
#' @param vectorB the second vector
#'
#' @return a list with 4 vectors, "only_in_vectorA" unique to vectorA, "only_in_vectorB" unique to vectorB, "shared_samples" shared between both vectorA and vectorB, "all" all values by combinng vectorA and vectorB 
set_decompose <- function(vectorA, vectorB){
  ret = list()
  # Find unique and shared samples
  ret[["only_in_vectorA"]] <- setdiff(vectorA, vectorB)  # Samples only in vectorA
  ret[["only_in_vectorB"]] <- setdiff(vectorB, vectorA)  # Samples only in vectorB
  ret[["shared_samples"]] <- intersect(vectorA, vectorB) # Samples shared between vectorA and vectorB
  ret[["all"]] <- union(vectorA, vectorB) # All samples between vectorA and vectorB
  return(ret)
}

#' Find missing columns from a tibble 
#'
#' @param tib the tibble to check
#' @param columns the columns to check for
#'
#' @return returns any missing columns 
get_missing_cols <- function(tib, columns) {
  setdiff(columns, colnames(tib))
}


#' Read and check microhaplotypes allele frequencies 
#'
#' @param mhaps_slaf_fnp 
#'
#' @returns the allele frequencies per microhaplotype 
#' @export
#'
process_input_mhaps_slaf <- function(mhaps_slaf_fnp){
  mhaps_slaf = readr::read_tsv(mhaps_slaf_fnp)
  required_mhaps_columns = c("target_name", "seq", "freq", "sample_total")
  input = readr::read_tsv(mhaps_slaf_fnp)
  missing_cols = get_missing_cols(input, required_mhaps_columns)
  if(length(missing_cols) > 0){
    stop(paste0("missing the following columns: ", paste0(required_mhaps_columns, collapse = ", "), " from ", mhaps_slaf_fnp) )
  }
  # validate columns mhaps_slaf_fnp
  mhaps_slaf_rules <- validate::validator(
    is.character(target_name), 
    is.character(seq),
    is.numeric(freq),
    is.numeric(sample_total), 
    ! is.na(target_name), 
    ! is.na(seq), 
    ! is.na(freq), 
    ! is.na(sample_total)
  )
  mhaps_slaf_fails <- validate::confront(input, mhaps_slaf_rules, raise = "all") %>%
    validate::summary() %>%
    dplyr::filter(fails > 0)
  
  if (nrow(mhaps_slaf_fails) > 0) {
    warns = paste0(
      "Input ", mhaps_slaf_fnp, " failed one or more validation checks: ", 
      str_c(mhaps_slaf_fails$expression, collapse = "\n")
    )
  }
  return(input)
}


#' Read and check translated loci of interest associated with the microhaplotype sequence  
#'
#' @param loci_of_interest_per_microhaps_fnp 
#'
#' @returns translated loci of interest  
#' @export
#'
process_input_loci_of_interest_per_microhaps <- function(loci_of_interest_per_microhaps_fnp){
  loci_of_interest_per_microhaps = readr::read_tsv(loci_of_interest_per_microhaps_fnp)
  required_mhaps_columns = c("target_name", "gene_id", "aa_position", "seq", "aa")
  input = readr::read_tsv(loci_of_interest_per_microhaps_fnp)
  missing_cols = get_missing_cols(input, required_mhaps_columns)
  if(length(missing_cols) > 0){
    stop(paste0("missing the following columns: ", paste0(required_mhaps_columns, collapse = ", "), " from ", loci_of_interest_per_microhaps_fnp) )
  }
  # validate columns loci_of_interest_per_microhaps_fnp
  loci_of_interest_per_microhaps_rules <- validate::validator(
    is.character(target_name), 
    is.character(gene_id),
    is.numeric(aa_position),
    is.character(seq), 
    is.character(aa),
    ! is.na(target_name), 
    ! is.na(gene_id), 
    ! is.na(aa_position), 
    ! is.na(seq), 
    ! is.na(aa)
  )
  loci_of_interest_per_microhaps_fails <- validate::confront(input, loci_of_interest_per_microhaps_rules, raise = "all") %>%
    validate::summary() %>%
    dplyr::filter(fails > 0)
  
  if (nrow(loci_of_interest_per_microhaps_fails) > 0) {
    warns = paste0(
      "Input ", loci_of_interest_per_microhaps_fnp, " failed one or more validation checks: ", 
      str_c(loci_of_interest_per_microhaps_fails$expression, collapse = "\n")
    )
  }
  return(input)
}

# Parse arguments ------------------------------------------------------
opts <- list(
  make_option(
    "--mhaps_slaf_fnp", 
    help = str_c(
      "TSV containing the columns: target_name, seq, freq, sample_total. The target_name and seq columns should match up with the columns in loci_of_interest_per_microhaps_fnp"
    )
  ),
  make_option(
    "--loci_of_interest_per_microhaps_fnp", 
    help = str_c(
      "TSV containing the columns: target_name, seq, gene_id, aa_position, aa. The target_name and seq columns should match up with the columns in mhaps_slaf_fnp"
    )
  ), 
  make_option(
    "--slaf_output", 
    help = str_c(
      "the output for the single locus allele frequency, will collapse frequencies across overlapping targets"
    )
  ), 
  make_option(
    "--per_target_slaf_output", 
    help = str_c(
      "optional output for the single lcous allele frequency calculated per target"
    )
  )
)
# parse arguments
args <- parse_args(OptionParser(option_list = opts))

## check for required arguments
required_arguments = c("mhaps_slaf_fnp", "loci_of_interest_per_microhaps_fnp", "slaf_output")
# checkOptparseRequiredArgsThrow(args, required_arguments)

if(interactive()){
  args$mhaps_slaf_fnp = "../../data/example_mhaps_slaf.tsv" 
  args$loci_of_interest_per_microhaps_fnp = "../../data/example_loci_of_interest_for_target_for_microhap.tsv" 
  args$slaf_output = "slaf.tsv" 
  args$per_target_slaf_output = "per_target_slaf.tsv" 
}

# read and check input 
mhaps_slaf = process_input_mhaps_slaf(args$mhaps_slaf_fnp)
translated_mhaps = process_input_loci_of_interest_per_microhaps(args$loci_of_interest_per_microhaps_fnp)

# join together tables 
# the translated mhaps may come from a full population and therefore there might be translated seqs that are missing from 
# the population frequencies and there might be hap frequencies that haven't been translated so will do inner join 
combined_tables = mhaps_slaf  %>% 
  inner_join(translated_mhaps, by = c("target_name", "seq"))

# calculate per target, renormalize freq in case input's freqs do not add up to 1 
translated_mhaps_slaf_per_target = combined_tables %>% 
  group_by(target_name, gene_id, aa_position, sample_total, aa) %>% 
  summarise(freq = sum(freq)) %>% 
  group_by(target_name, gene_id, aa_position, sample_total) %>% 
  mutate(total_freq = sum(freq)) %>% 
  mutate(freq = freq/total_freq) %>% 
  select(-total_freq) %>% 
  unite(variant, gene_id, aa_position, aa, sep = ":")

# calculate collapsed over target, will weight evenly between overlapping targets
# @todo consider weighting by the sample total per target to give sample weighted in case one target has very poor coverage 
# renormalize freq in case input's freqs do not add up to 1
# take the max sample total between targets to get the sample total 
translated_mhaps_slaf = combined_tables %>% 
  group_by(gene_id, aa_position, aa) %>% 
  summarise(freq = sum(freq), 
            sample_total = max(sample_total)) %>% 
  group_by(gene_id, aa_position, sample_total) %>% 
  mutate(total_freq = sum(freq)) %>% 
  mutate(freq = freq/total_freq) %>% 
  select(-total_freq) %>% 
  unite(variant, gene_id, aa_position, aa, sep = ":")

# write output 
write_tsv(translated_mhaps_slaf, args$slaf_output)

# optionally write output per target if output name given 
if(!is.null(args$per_target_slaf_output)){
  write_tsv(translated_mhaps_slaf_per_target, args$per_target_slaf_output)
}
