#!/usr/bin/env Rscript


packagesToLoad = c("tibble", "dplyr", "stringr", "readr", "optparse", "Biostrings")

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
    "--ref_bed", 
    help = str_c(
      "a bed file containing the reference location of the ref_seq, no column names but the first 6 columns should be chrom, start, end, target_id, length, strand"
    )
  ), 
  make_option(
    "--fasta", 
    help = str_c(
      "a fasta file with the ref sequences, the names of the records should match up with the target_id of the --ref_bed file"
    )
  ), 
  make_option(
    "--out", 
    help = str_c(
      "the out file to write to"
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

required_arguments = c("ref_bed", "out")




# parse arguments
arg <- parse_args(OptionParser(option_list = opts))
## check for required arguments
checkOptparseRequiredArgsThrow(arg, required_arguments)

# read in input data and gather warnings about columns and data formatting 
warnings = c()

# read in the panel reference information 
ref_bed = readr::read_tsv(arg$ref_bed, col_names = T)
warnings = c(warnings, genWarningsMissCols(ref_bed, c("#chrom", "start", "end", "target_id", "length", "strand"), arg$ref_bed))

if(file.exists(arg$out) & !arg$overwrite){
  warnings = c(warnings, "file ", arg$out, " already exist, use --overwrite to over write it") 
}



dna <- Biostrings::readDNAStringSet(arg$fasta)


# check to see if values between dataets are similar 
ref_allele_decomp = set_decompose(ref_bed$target_id, unique(allele_table$target_id))

if(length(ref_allele_decomp$only_in_vectorA) > 0){
  warnings = c(warnings, paste0("the following loci were missing from the fasta file ", arg$fasta, " but are in ", arg$ref_bed,
                                "\n", paste0(ref_allele_decomp$only_in_vectorA, collapse = ",")
  ) 
  )
}


if(length(warnings) > 0){
  stop(paste0("\n", paste0(warnings, collapse = "\n")) )
}











