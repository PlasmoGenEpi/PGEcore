#!/usr/bin/env Rscript


packagesToLoad = c("tibble", "dplyr", "stringr", "readr", "optparse")

loaded = suppressMessages(lapply(packagesToLoad, require, character.only = TRUE))

options(readr.show_col_types = FALSE)
options(dplyr.summarise.inform = FALSE)

 
#' Check for required arguments, and report which are missing 
#'
#' @param arg the parsed arguments from arg parse
#' @param required_args the required arguments
#'
#' @returns returns void if all required arguments
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
#' @returns a list with 4 vectors, "only_in_vectorA" unique to vectorA, "only_in_vectorB" unique to vectorB, "shared_samples" shared between both vectorA and vectorB, "all" all values by combinng vectorA and vectorB 
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
#' @returns returns any missing columns 
returnMissingColumns <-function(tib, columns){
  setdiff(columns, colnames(tib))
}


#' Add a warning about Missing columns 
#'
#' @param tib the tibble to check
#' @param cols the columns to check for
#' @param fnp the file name path from which the tibble was read from
#'
#' @returns adds a warning message about any missing columns 
genWarningsMissCols <- function(tib, cols, fnp){
  col_miss = returnMissingColumns(tib, cols)
  if(length(col_miss)> 0){
    return(paste0("missing the following columns: ", paste0(col_miss, collapse = ","), " from ", fnp))
  }
  return (character(0))
}


# Parse arguments ------------------------------------------------------
opts <- list(
  make_option(
    "--amino_acid_calls", 
    help = str_c(
      "table with amino acid calls"
    )
  ),
  make_option(
    "--out", 
    help = str_c(
      "the out file to write to"
    )
  ), 
  make_option(
    "--out_nonbiallelic", 
    help = str_c(
      "only process these targets"
    )
  ), 
  make_option(
    "--overwrite", 
    help = str_c(
      "overwrite the output if it already exists"
    ), 
    action="store_true", 
    default=FALSE,
    type = "logical"
  )
)

required_arguments = c("amino_acid_calls", "out")




# parse arguments
arg <- parse_args(OptionParser(option_list = opts))
## check for required arguments
checkOptparseRequiredArgsThrow(arg, required_arguments)

amino_acid_calls = readr::read_tsv(arg$amino_acid_calls)

warnings = genWarningsMissCols(amino_acid_calls, c("gene_id", "aa_position", "ref_aa", "aa"))

if(file.exists(arg$out) & !arg$overwrite){
  warnings = c(warnings, "file ", arg$out, " already exist, use --overwrite to over write it") 
}

if("out_nonbiallelic" %in% names(arg) && file.exists(arg$out_nonbiallelic) && !arg$overwrite){
  warnings = c(warnings, "file ", arg$out_nonbiallelic, " already exist, use --overwrite to over write it") 
}

if(length(warnings) > 0){
  stop(paste0("\n", paste0(warnings, collapse = "\n")) )
}

amino_acid_calls = amino_acid_calls |> 
  group_by(gene_id, aa_position, ref_aa) |> 
  mutate(allele_calls = n_distinct(aa))

amino_acid_calls_biallelic = amino_acid_calls |> 
  filter(allele_calls <= 2)

amino_acid_calls_more_than_biallelic = amino_acid_calls |> 
  filter(allele_calls > 2)

write_tsv(amino_acid_calls_biallelic, arg$out)
if("out_nonbiallelic"  %in% names(arg)){
  write_tsv(amino_acid_calls_more_than_biallelic, arg$out_nonbiallelic)
}

