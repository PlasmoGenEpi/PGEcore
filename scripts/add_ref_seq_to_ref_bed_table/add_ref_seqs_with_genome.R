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


#' Gen warns if ref_bed columns are not the expected types 
#'
#' @param ref_bed 
#'
#' @returns warnings about types and missingness 
gen_warns_on_validate_ref_bed_cols <- function(ref_bed){
  # validate columns
  rules <- validate::validator(
    is.character(`#chrom`), 
    is.integer(start),
    is.integer(end),
    is.character(target_id), 
    is.integer(length), 
    is.character(strand), 
    ! is.na(`#chrom`), 
    ! is.na(start), 
    ! is.na(end), 
    ! is.na(target_id), 
    ! is.na(length), 
    ! is.na(strand)
  )
  fails <- validate::confront(ref_bed, rules, raise = "all") %>%
    validate::summary() %>%
    dplyr::filter(fails > 0)
  warns = c()
  if (nrow(fails) > 0) {
    warns = paste0(
      "Input input_data failed one or more validation checks: ", 
      str_c(fails$expression, collapse = "\n")
    )
  }
  return(warns)
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
    "--genome", 
    help = str_c(
      "a genome file to extract the ref_seq"
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

#' Run adding reference seqs to a file by using the genome described in the reference bed file 
#'
#' @returns true if runs all the way through 
run_add_ref_seqs_with_genome_to_ref_bed <-function(){
  required_arguments = c("ref_bed", "genome", "out")
  
  # parse arguments
  arg <- parse_args(OptionParser(option_list = opts))
  ## check for required arguments
  checkOptparseRequiredArgsThrow(arg, required_arguments)
  
  # read in input data and gather warnings about columns and data formatting 
  warnings = c()
  
  # read in the panel reference information 
  ref_bed = readr::read_tsv(arg$ref_bed, col_names = T)
  
  warnings = c(warnings, genWarningsMissCols(ref_bed, c("#chrom", "start", "end", "target_id", "length", "strand"), arg$ref_bed))
  
  # check output options 
  if(file.exists(arg$out) & !arg$overwrite){
    warnings = c(warnings, "file ", arg$out, " already exist, use --overwrite to over write it") 
  }
  if(length(warnings) > 0){
    stop(paste0("\n", paste0(warnings, collapse = "\n")) )
  }
  warnings = c()

  # check columns types 
  warnings = c(warnings, gen_warns_on_validate_ref_bed_cols(ref_bed))
  
  # check if multiple target_id loaded 
  ref_bed_name_sum_multi = ref_bed |> 
    group_by(target_id) |> 
    count() |> 
    filter(n > 1)
  if(nrow(ref_bed_name_sum_multi) > 0){
    warnings = c(warnings, paste0("found multi names for target_id in ", arg$ref_bed, " found the following multiple times: ", paste0(ref_bed_name_sum_multi$target_id, collapse = ",")) ) 
  }
  
  if(length(warnings) > 0){
    stop(paste0("\n", paste0(warnings, collapse = "\n")) )
  }
  
  # read in gnome 
  loaded_genome = Biostrings::readDNAStringSet("/tank/data/genomes/plasmodium/genomes/pf/genomes/Pf3D7.fasta")
  # remove whitespace and beyond in names 
  names(loaded_genome) = gsub(" .*", "", names(loaded_genome))
  
  ref_bed$ref_seq = ""
  
  for(row in 1:nrow(ref_bed)){
    if(ref_bed$`#chrom`[row] %in% names(loaded_genome)){
      ref_seq = Biostrings::subseq(loaded_genome[ref_bed$`#chrom`[row]], ref_bed$start[row] + 1, width = ref_bed$length[row])
      if('-' == ref_bed$`#chrom`[row]){
        ref_seq = Biostrings::reverseComplement(ref_seq)
      }
      ref_bed$ref_seq[row] = ref_seq
    } else { 
      stop(paste0(ref_bed$`#chrom`[row], "not in ", arg$genome, " options: ", paste0(names(loaded_genome),collapse = ",")))
    }
  }
  
  write_tsv(ref_bed, arg$out)
  return (T)
}

run_add_ref_seqs_with_genome_to_ref_bed_status = run_add_ref_seqs_with_genome_to_ref_bed()
