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
#' @return returns void if all required arguments
checkOptparseRequiredArgsThrow <- function(arg, required_args){
  missing <- setdiff(required_args, names(arg))
  if(length(missing) > 0){
    missing = paste0("--", missing)
    stop(paste0("mssing the following arguments: ", paste0(missing, collapse = ", ")))
  }
}


#' Find missing columns from a tibble 
#'
#' @param tib the tibble to check
#' @param columns the columns to check for
#'
#' @return returns any missing columns 
returnMissingColumns <-function(tib, columns){
  setdiff(columns, colnames(tib))
}


#' Add a warning about Missing columns 
#'
#' @param tib the tibble to check
#' @param cols the columns to check for
#' @param fnp the file name path from which the tibble was read from
#'
#' @return adds a warning message about any missing columns 
genWarningsMissCols <- function(tib, cols, fnp){
  col_miss = returnMissingColumns(tib, cols)
  if(length(col_miss)> 0){
    return(paste0("missing the following columns: ", paste0(col_miss, collapse = ","), " from ", fnp))
  }
  return (character(0))
}



#' Check for sub selecting arguments of target_ids and specimen_ids
#'
#' @param parsed_args 
#'
#' @return a list with two named vectors, select_target_ids and select_specimen_ids
check_subselecting_args <- function(parsed_args){
  select_target_ids = c()
  if("select_target_ids" %in% names(parsed_args)){
    if(file.exists(parsed_args$select_target_ids)){
      select_target_ids = readr::read_tsv(parsed_args$select_target_ids, col_names = F)$X1
    } else { 
      select_target_ids = unlist(strsplit(parsed_args$select_target_ids, split = ","))
    }
  }
  
  select_specimen_ids = c()
  if("select_specimen_ids" %in% names(parsed_args)){
    if(file.exists(parsed_args$select_specimen_ids)){
      select_specimen_ids = readr::read_tsv(parsed_args$select_specimen_ids, col_names = F)$X1
    } else { 
      select_specimen_ids = unlist(strsplit(parsed_args$select_specimen_ids, split = ","))
    }
  }
  list(select_target_ids = select_target_ids, 
       select_specimen_ids = select_specimen_ids)
}


#' Fitler an snp table for select target_ids and specimen_ids 
#'
#' @param snp_data the snp data 
#' @param opt_sub_sels the optional subselections, a list with two named vectors: select_target_ids and select_specimen_ids
#'
#' @return a filtered snp_data table 
filter_snp_table_for_optional_subselecting <- function (snp_data, opt_sub_sels){
  if(length(opt_sub_sels$select_specimen_ids) > 0 ){
    snp_data = snp_data |> 
      filter(specimen_id %in% opt_sub_sels$select_specimen_ids)
  }
  
  if(length(opt_sub_sels$select_target_ids) > 0 ){
    snp_data = snp_data |> 
      filter(target_id %in% opt_sub_sels$select_target_ids)
  }
  return(snp_data)
}


#' generate warnings for sub selecting for selections that don't exist 
#'
#' @param snp_data the snp data 
#' @param opt_sub_sels the optional subselections, a list with two named vectors: select_target_ids and select_specimen_ids 
#' @param snp_table_fnp the table the snp data was read from 
#'
#' @return warnings about missing sub-selections
check_warnings_for_subselecting_snp_table <- function(snp_data, opt_sub_sels, snp_table_fnp){
  warns = c()
  if(length(opt_sub_sels$select_specimen_ids) > 0 ){
    missing_sel_specs = setdiff(opt_sub_sels$select_specimen_ids, unique(snp_data$specimen_id))
    if(length(missing_sel_specs) > 0 ){
      warns = c(warns, paste0("supplied --select_specimen_ids but the following specimen_ids are missing from ", snp_table_fnp, "\n", 
                              paste0(missing_sel_specs, collapse = ",")))
    }
  }
  
  if(length(opt_sub_sels$select_target_ids) > 0 ){
    missing_sel_tars = setdiff(opt_sub_sels$select_target_ids, unique(snp_data$target_id))
    if(length(missing_sel_tars) > 0 ){
      warns = c(warns, paste0("supplied --select_target_ids but the following target_ids are missing from ", snp_table_fnp, "\n", 
                              paste0(missing_sel_tars, collapse = ",")))
    }
  }
  return (warns)
}

# Parse arguments ------------------------------------------------------
opts <- list(
  make_option(
    "--snp_table", 
    help = str_c(
      "TSV containing at least the columns: specimen_id, target_id, chrom, pos, snp_name, ref_base, seq_base, read_count, is_biallelic"
    )
  ),
  make_option(
    "--output_fnp", 
    help = str_c(
      "the output filename path for the output filtered file"
    )
  ),
  make_option(
    "--mindist_between_snps", 
    type = "integer",
    default = 10000,
    help = str_c(
      "the minimum distance between SNPs to consider them 'independent'. Default: %default"
    )
  ), 
  make_option(
    "--select_target_ids", 
    help = str_c(
      "only process these targets"
    )
  ), 
  make_option(
    "--select_specimen_ids", 
    help = str_c(
      "only process these samples"
    )
  ), 
  make_option(
    "--overwrite", 
    help = str_c(
      "overwrite the output fnp if it already exists. Default: %default"
    ), 
    action="store_true", 
    default=FALSE,
    type = "logical"
  ), 
  make_option(
    "--only_biallelic", 
    help = str_c(
      "filter to only biallelic SNPs too. Default: %default"
    ), 
    action="store_true", 
    default=FALSE,
    type = "logical"
  ), 
  make_option(
    "--only_informative", 
    help = str_c(
      "filter off SNPs with no variation. Default: %default"
    ), 
    action="store_true", 
    default=FALSE,
    type = "logical"
  )
)


#' validate the input data columns 
#'
#' @param snp_table the snp data table 
#'
#' @return the warnings about the input 
validate_columns_types <-function(snp_table){
  warns = c()
  # validate columns snp_table
  snp_table_rules <- validate::validator(
    is.character(specimen_id),
    is.numeric(read_count),
    is.character(target_id), 
    is.character(ref_base),
    is.character(seq_base), 
    is.character(chrom),
    is.character(snp_name),
    is.numeric(pos),
    is.logical(is_biallelic),
    ! is.na(specimen_id), 
    ! is.na(read_count), 
    ! is.na(target_id), 
    ! is.na(ref_base),
    ! is.na(seq_base), 
    ! is.na(chrom),
    ! is.na(snp_name),
    ! is.na(pos),
    ! is.na(is_biallelic)
  )
  snp_table_fails <- validate::confront(snp_table, snp_table_rules, raise = "all") |>
    validate::summary() |>
    dplyr::filter(fails > 0)
  if (nrow(snp_table_fails) > 0) {
    warns = paste0(
      "Input snp_table failed one or more validation checks: ", 
      str_c(snp_table_fails$expression, collapse = "\n")
    )
  }
  return(warns)
}

#' run filtering SNPs to the highest diversity SNP per microhaplotype target
#'
#' @return true if runs all the way through
run_filter_to_highest_diversity_snp_call <-function(){
  required_arguments = c("snp_table", "output_fnp")
  # parse arguments
  arg <- parse_args(OptionParser(option_list = opts))
  
  ## check for required arguments
  checkOptparseRequiredArgsThrow(arg, required_arguments)
  
  # check if output exists and if it does throw if not overwriting 
  if(file.exists(arg$output_fnp) & !arg$overwrite){
    stop(paste0("file ", arg$output_fnp, " exits, use --overwrite to over write it"))
  }
  
  # process if only processing specific targets and specimens
  optional_sub_selections = check_subselecting_args(arg)
  
  # read in input data and gather warnings about columns and data formatting 
  warnings = c()

  # read in snp table for the microhaplotype data 
  snp_table = readr::read_tsv(arg$snp_table)
  
  # sanity check the input
  warnings = c(warnings, genWarningsMissCols(snp_table, c("specimen_id", "target_id", "chrom", "pos", "snp_name", "ref_base", "seq_base", "read_count", "is_biallelic"), arg$snp_table))
  if(length(warnings) > 0){
    stop(paste0("\n", paste0(warnings, collapse = "\n")) )
  }
  warnings = c()
  warnings = c(warnings, check_warnings_for_subselecting_snp_table(snp_table, optional_sub_selections, arg$snp_table))
  if(length(warnings) > 0){
    stop(paste0("\n", paste0(warnings, collapse = "\n")) )
  }
  
  # filter snp table for optional sub-selecting
  snp_table = filter_snp_table_for_optional_subselecting(snp_table, optional_sub_selections)
  
  # filter to biallelic SNPs too if set 
  if(arg$only_biallelic){
    snp_table = snp_table |>
      filter(is_biallelic)
  }
  
  # for counts and he calc to work, SNPs need to be collapsed already 
  snp_table_counts_per_specimen_id = snp_table |> 
    group_by(specimen_id, chrom, pos, snp_name, ref_base, seq_base) |> 
    mutate(n = n())
  
  warnings = c()
  snp_table_counts_per_specimen_id_multi = snp_table_counts_per_specimen_id |> 
    filter(n > 1)
  if(nrow(snp_table_counts_per_specimen_id_multi) > 0){
    warnings = c(warnings, paste0("the following snps were found to have multiple calls per specimen_id for the same seq_base, make sure calls are collapsed", 
                                  ":\n", 
                                  paste0(unique(snp_table_counts_per_specimen_id_multi$snp_name), collapse = ",")))
  }
  if(length(warnings) > 0){
    stop(paste0("\n", paste0(warnings, collapse = "\n")) )
  }
  
  # get expected heterozygosity (he) of SNPs 
  snp_table_freqs = snp_table |> 
    group_by(chrom, pos, snp_name, ref_base, seq_base) |> 
    summarise(allele_count = n()) |> 
    group_by(chrom, pos, snp_name, ref_base) |> 
    mutate(allele_total = sum(allele_count)) |> 
    mutate(allele_freq = allele_count/allele_total)
  
  snp_table_he = snp_table_freqs |> 
    group_by(chrom, pos, snp_name, ref_base) |> 
    summarise(he = 1 - sum(allele_freq^2)) |> 
    arrange(desc(he))
  
  # filtering table by getting max he and then taking the next best he that doesn't come within the min distance
  snp_table_he$keep = T
  for(row in 2:nrow(snp_table_he)){
    for(filt_row in 1:(row -1)){
      if(snp_table_he$chrom[filt_row] == snp_table_he$chrom[row] & 
         abs(snp_table_he$pos[filt_row] - snp_table_he$pos[row]) < arg$mindist_between_snps){
        snp_table_he$keep[row] = F
        break
      }
    }
  }
  snp_table_he_filt = snp_table_he |> 
    filter(keep) |> 
    mutate(snp_id = paste0(chrom, "-", pos))
  
  if(arg$only_informative){
    # if turned on, filter only to informative SNPs e.g. those with variation
    snp_table_he_filt = snp_table_he_filt |> 
      filter(he > 0)
  }
  # join back and get the max per target_id 
  snp_table = snp_table |> 
    left_join(snp_table_he |> 
                select(-keep), by = c("chrom", "pos", "snp_name", "ref_base"))
  
  # filter to the max he with no ties (first come first out slicing)
  # to get the best snp for each target
  snp_table_out = snp_table |> 
    filter(paste0(chrom, "-", pos)  %in%  snp_table_he_filt$snp_id)

  # writing output results 
  write_tsv(snp_table_out, arg$output_fnp)
  
  return(T)
}

run_status = run_filter_to_highest_diversity_snp_call()
