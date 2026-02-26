#!/usr/bin/env Rscript 

# quietly load packages 
packagesToLoad = c("dplyr", "optparse", "purrr", "readr", "stringr", "tidyr", "validate")

loaded = suppressMessages(lapply(packagesToLoad, require, character.only = TRUE))

options(readr.show_col_types = FALSE)
options(dplyr.summarise.inform = FALSE)

# function to get not in (makes it easier than doing !(test %in% test_set) )
`%!in%` <- Negate(`%in%`)

# Parse arguments ------------------------------------------------------
opts <- list(
  make_option(
    "--aa_table", 
    type = "character",
    help = str_c(
      "Path to a TSV file containing amino acid calls, with the columns: specimen_id, ", 
      "gene, pos, read_count, aa"
    )
  ), 
  make_option(
    "--loci_groups_input", 
    type = "character",
    help = str_c(
      "Path to a TSV file containing loci group definitions, with the ", 
      "columns: group_id, gene_id, aa_position"
    )
  ), 
  make_option(
    "--output_path",
    type = "character",
    help = str_c(
      "Path to write an output TSV file containing prevalence and frequency estimates ",
      "for all variants in the data"
    )
  ), 
  make_option(
    "--recalc_single_locus_output_path",
    type = "character",
    help = str_c(
      "An optional output for the re-calculated single locus frequences/prevalences from the multilocus, ",
      "can be helpful for checking if the multilocus is capturing well all available variants"
    )
  ),
  make_option(
    c("--method"),
    help = stringr::str_c(
      "Method to use for estimating allele frequency. Options are: ",
      "wsaf_prop (frequencies are weighted by their within-sample-allele-frequency, no allele_counts will be exported in this case), presence_absence (simply count by how many times each allele appears). Default: %default"
    ),
    type = "character",
    default = "wsaf_prop",
    callback = function(opt, flag_string, value, parser, ...) {
      if (!value %in% c("wsaf_prop", "presence_absence")) {
        stop(stringr::str_c(value, " is not a valid method"))
      }
      value
    }
  ),
  make_option(
    "--wsaf_cut_off",
    type = "double",
    default = 0.70,
    help = str_c(
      "To help determine multi-locus variants, we filter off minor variants to determine ", 
      "a \"dominant\" strain at each loci that will then be reasonable to assume at least ", 
      "one haplotype with all variants above that cut off exists"
    )
  )
)

# Set up options
#' Check for required arguments, and report which are missing 
#'
#' @param parser the parser created from optparse
#' @param arg the parsed arguments from optparse
#' @param required_args the required arguments (without the --)
#'
#' @return returns void if all required arguments
checkOptparseRequiredArgsThrow <- function(parser, arg, required_args){
  missing <- setdiff(required_args, names(arg))
  if(length(missing) > 0){
    missing = paste0("--", missing)
    print_help(parser)
    stop(paste0("mssing the following arguments: ", paste0(missing, collapse = ", ")))
  }
}

#' Read amino acid calls table into a tibble
#'
#' @description
#' Takes the path of a TSV file of amino acid calls. Reads these in, 
#' computes the number of amino acids at each position, and returns a 
#' tibble with this information.
#'
#' @param aa_table Path to a TSV file with amino acid calls. It should have
#'   columns for specimen_id, gene, pos, read_count, aa.
#' 
#' @import dplyr
#' 
#' @return Tibble of amino acid calls with specimen_id, gene, pos, 
#'   read_count, aa, and n_aa columns.
create_aa_table_input <- function(aa_table) {
  
  # Check input arguments
  stopifnot(is.character(aa_table))
  
  # read in amino acid calls and validate columns
  df_aa <- read.table(aa_table, header = TRUE)
  rules <- validate::validator(
    is.character(specimen_id), 
    is.character(gene), 
    is.character(gene_id),
    is.integer(aa_position), 
    is.integer(read_count), 
    is.character(aa), 
    ! is.na(specimen_id), 
    ! is.na(gene), 
    ! is.na(gene_id), 
    ! is.na(aa_position), 
    ! is.na(read_count), 
    ! is.na(aa)
  )
  fails <- validate::confront(df_aa, rules, raise = "all") %>%
    validate::summary() %>%
    dplyr::filter(fails > 0)
  if (nrow(fails) > 0) {
    stop(
      "Input input_data failed one or more validation checks: ", 
      str_c(fails$expression, collapse = "\n"), 
      call. = FALSE
    )
  }
  
  # get number of amino acids at each locus
  df_aa <- df_aa |>
    group_by(specimen_id, gene_id, aa_position) |>
    mutate(n_aa = n()) |>
    ungroup()

  return(df_aa)
}

#' Read in loci groups table
#'
#' This function takes the path to a table of loci groups for 
#' multilocus allele frequency and prevalence calculations and reads it 
#' into a tibble.
#'
#' @param loci_groups_path Path to loci groups TSV. It should have 
#'   columns for group_id, gene_id, and aa_position.
#'
#' @import dplyr
#'
#' @return Tibble of loci groups, with columns for group_id, gene_id, 
#'   and aa_position.
create_loci_group_input <- function(loci_groups_path) {

  # Check input arguments
  stopifnot(is.character(loci_groups_path))
  
  # Read and validate table
  loci_groups <- read_tsv(
    loci_groups_path, 
    col_types = cols(.default = col_character(), aa_position = col_integer()), 
    progress = FALSE
  )
  rules <- validate::validator(
    is.character(group_id), 
    is.character(gene_id), 
    is.integer(aa_position), 
    ! is.na(group_id), 
    ! is.na(gene_id), 
    ! is.na(aa_position)
  )
  fails <- validate::confront(loci_groups, rules, raise = "all") %>%
    validate::summary() %>%
    dplyr::filter(fails > 0)
  if (nrow(fails) > 0) {
    stop(
      "Input input_data failed one or more validation checks: ", 
      str_c(fails$expression, collapse = "\n"), 
      call. = FALSE
    )
  }

  return(loci_groups)

}


#' Generate the single locus allele frequencies and prevalence from a table of multilocus calls weighted by wsaf 
#'
#' @param multilocus_calls A table of 3 columns, group_id, specimen_id, and variant (in variantstring format)
#'
#' @returns summarized per group, single locus freq/prev recalculated from the multilocus variant calls 
#' @export
#'
#' @examples
generate_single_locus_prev_freq_from_multilocus_groups_wsaf_prop <-function(multilocus_calls){
  # validate input
  # group_id, specimen_id, variant 
  rules <- validate::validator(
    is.character(group_id), 
    is.character(specimen_id), 
    is.character(variant),
    is.numeric(wsaf),
    ! is.na(group_id), 
    ! is.na(specimen_id), 
    ! is.na(variant), 
    ! is.na(wsaf)
  )
  fails <- validate::confront(multilocus_calls, rules, raise = "all") %>%
    validate::summary() %>%
    dplyr::filter(fails > 0)
  if (nrow(fails) > 0) {
    stop(
      "Input multilocus_calls for function generate_single_locus_prev_freq_from_multilocus_groups failed one or more validation checks: ", 
      str_c(fails$expression, collapse = "\n"), 
      call. = FALSE
    )
  }
  # deconstruct the variantstring format into single locus 
  multilocus_calls_mod = multilocus_calls %>% 
    group_by(group_id) %>% 
    mutate(variant_split = strsplit(variant, split = ";")) %>% 
    unnest(variant_split) %>% 
    separate(variant_split, into = c("gene_id", "aa_position", "aa"), sep = ":")%>% 
    mutate(aa_position = strsplit(aa_position, split = "_"), 
           aa = strsplit(aa, split = "_")) %>% 
    unnest(c(aa_position, aa)) %>% 
    mutate(aa_position = as.numeric(aa_position)) %>% 
    arrange(group_id, gene_id, aa_position, aa)
  
  # now re-calculate allele frequency and prev 
  multilocus_calls_mod_slaf = multilocus_calls_mod %>% 
    group_by(group_id, gene_id, aa_position) %>% 
    mutate(sample_total = n_distinct(specimen_id), 
           wsaf_total = sum(wsaf)) %>% 
    group_by(group_id, gene_id, aa_position, aa, sample_total, wsaf_total) %>% 
    summarise(wsaf_sum = sum(wsaf), 
              sample_count = n_distinct(specimen_id)) %>% 
    mutate(freq = wsaf_sum/wsaf_total, 
           prev = sample_count/sample_total) %>% 
    unite(variant,gene_id, aa_position, aa, sep = ":") %>% 
    select(-wsaf_sum,-wsaf_total)
  
  return (multilocus_calls_mod_slaf)
}

#' Generate the single locus allele frequencies and prevalence from a table of multilocus calls based on presence absence
#'
#' @param multilocus_calls A table of 3 columns, group_id, specimen_id, and variant (in variantstring format)
#'
#' @returns summarized per group, single locus freq/prev recalculated from the multilocus variant calls 
#' @export
#'
#' @examples
generate_single_locus_prev_freq_from_multilocus_groups_presence_absence <-function(multilocus_calls){
  # validate input
  # group_id, specimen_id, variant 
  rules <- validate::validator(
    is.character(group_id), 
    is.character(specimen_id), 
    is.character(variant), 
    ! is.na(group_id), 
    ! is.na(specimen_id), 
    ! is.na(variant)
  )
  fails <- validate::confront(multilocus_calls, rules, raise = "all") %>%
    validate::summary() %>%
    dplyr::filter(fails > 0)
  if (nrow(fails) > 0) {
    stop(
      "Input multilocus_calls for function generate_single_locus_prev_freq_from_multilocus_groups failed one or more validation checks: ", 
      str_c(fails$expression, collapse = "\n"), 
      call. = FALSE
    )
  }
  # deconstruct the variantstring format into single locus 
  multilocus_calls_mod = multilocus_calls %>% 
    group_by(group_id) %>% 
    mutate(variant_split = strsplit(variant, split = ";")) %>% 
    unnest(variant_split) %>% 
    separate(variant_split, into = c("gene_id", "aa_position", "aa"), sep = ":")%>% 
    mutate(aa_position = strsplit(aa_position, split = "_"), 
           aa = strsplit(aa, split = "_")) %>% 
    unnest(c(aa_position, aa)) %>% 
    mutate(aa_position = as.numeric(aa_position)) %>% 
    arrange(group_id, gene_id, aa_position, aa)
  
  # now re-calculate allele frequency and prev 
  multilocus_calls_mod_slaf = multilocus_calls_mod %>% 
    group_by(group_id, gene_id, aa_position) %>% 
    mutate(sample_total = n_distinct(specimen_id), 
           allele_total = n()) %>% 
    group_by(group_id, gene_id, aa_position, aa, sample_total, allele_total) %>% 
    summarise(allele_count = n(), 
              sample_count = n_distinct(specimen_id)) %>% 
    mutate(freq = allele_count/allele_total, 
           prev = sample_count/sample_total) %>% 
    unite(variant,gene_id, aa_position, aa, sep = ":")
  
  return (multilocus_calls_mod_slaf)
}


#' Calculate the allele frequency and prevalence of multilocus allele weighted by their within sample frequency (wsaf)
#'
#' @param multilocus_calls - a table with columns specimen_id, variant, wsaf 
#'
#' @returns a tibble with freq and prev calculated 
#' @export
#'
calculate_multilocus_af_prev_wsaf_prop<-function(multilocus_calls){
  # validate input
  # specimen_id variant wsaf
  rules <- validate::validator(
    is.character(specimen_id), 
    is.character(variant), 
    is.numeric(wsaf),
    ! is.na(specimen_id), 
    ! is.na(variant), 
    ! is.na(wsaf)
  )
  fails <- validate::confront(multilocus_calls, rules, raise = "all") %>%
    validate::summary() %>%
    dplyr::filter(fails > 0)
  if (nrow(fails) > 0) {
    stop(
      "Input multilocus_calls for function calculate_multilocus_af_prev_presence_absence failed one or more validation checks: ", 
      str_c(fails$expression, collapse = "\n"), 
      call. = FALSE
    )
  }
  multilocus_calls_prev_freq = multilocus_calls %>% 
    ungroup() %>% 
    mutate(sample_total = n_distinct(specimen_id), 
           wsaf_total = sum(wsaf)) %>% 
    group_by(variant, sample_total, wsaf_total) %>% 
    summarise(wsaf_sum = sum(wsaf), 
              sample_count = n_distinct(specimen_id)) %>% 
    mutate(prev = sample_count/sample_total, 
           freq = wsaf_sum/wsaf_total) %>% 
    ungroup() %>% 
    select(-wsaf_sum,-wsaf_total)
  return(multilocus_calls_prev_freq)
}

#' Calculate the allele frequency and prevalence of multilocus allele by their presence/absence 
#'
#' @param multilocus_calls - a table with columns specimen_id, variant 
#'
#' @returns a tibble with freq and prev calculated 
#' @export
#'
calculate_multilocus_af_prev_presence_absence<-function(multilocus_calls){
  # validate input
  # specimen_id variant 
  rules <- validate::validator(
    is.character(specimen_id), 
    is.character(variant), 
    ! is.na(specimen_id), 
    ! is.na(variant)
  )
  fails <- validate::confront(multilocus_calls, rules, raise = "all") %>%
    validate::summary() %>%
    dplyr::filter(fails > 0)
  if (nrow(fails) > 0) {
    stop(
      "Input multilocus_calls for function calculate_multilocus_af_prev_presence_absence failed one or more validation checks: ", 
      str_c(fails$expression, collapse = "\n"), 
      call. = FALSE
    )
  }
  multilocus_calls_prev_freq = multilocus_calls %>% 
    ungroup() %>% 
    mutate(sample_total = n_distinct(specimen_id), 
           allele_total = n()) %>% 
    group_by(variant, sample_total, allele_total) %>% 
    summarise(allele_count = n(), 
              sample_count = n_distinct(specimen_id)) %>% 
    mutate(prev = sample_count/sample_total, 
           freq = allele_count/allele_total) %>% 
    ungroup()
  return(multilocus_calls_prev_freq)
}







# RUN MODULE -----------------------------------------------------------


# parse arguments
parser <- OptionParser(option_list = opts)
args <- parse_args(parser)
# Arguments used for development, in order to properly test the below, working directory has to be set the directory of this script
if (interactive()) {
  args$aa_table <- "data/example2_amino_acid_calls.tsv"
  args$loci_groups_input <- "data/example_loci_groups.tsv"
  args$output_path <- "mlafp.tsv"
  args$recalc_single_locus_output_path <- "recalc_sl_from_ml.tsv"
  args$wsaf_cut_off = 0.70
}

# make sure the required options are being provided 
required_arguments = c("aa_table", "loci_groups_input", "output_path")
checkOptparseRequiredArgsThrow(parser, args, required_arguments)

# Read in amino acid data 
aa_table <- create_aa_table_input(args$aa_table)

# Read in the loci groups to determine 
loci_groups <- create_loci_group_input(args$loci_groups_input) %>% 
  group_by(group_id) %>% 
  mutate(loci_in_group = n_distinct(paste0(gene_id, "-", aa_position)))

# split the loci into the groups that will then be processed 
loci_groups_split = split(loci_groups, loci_groups$group_id)

# create empty tables to fill with the final results 
all_aa_table_group_filt_final_prev_freq = tibble()
all_aa_table_group_filt_final = tibble()

for(loci_group in names(loci_groups_split)){
  # filter to loci in group and filter to only samples that have calls for all loci in the group 
  # the inner join here will make it so only the loci of the group of interest will be selected 
  aa_table_group = aa_table %>% 
    inner_join(loci_groups_split[[loci_group]], by = c("gene_id", "aa_position")) %>% 
    group_by(specimen_id) %>% 
    mutate(loci_called = n_distinct(paste0(gene_id, "-", aa_position))) %>% 
    filter(loci_called == loci_in_group) %>% 
    ungroup()
  
  # as a first pass will, will first create all the multi-locus groups for when only one of the loci is variable 
  # e.g. if have 3 loci, 1st and 2nd loci are A and 3rd loci has calls of RS, then we know that A-A-R and A-A-S exist
  # albeit this relies on reliable genotyping of all loci, if the 1st and 2nd loci don't have good depth and didn't capture minor 
  # alleles or the 3rd loci is a genotyping error than high chance of creating false haplotypes 
  aa_table_group_only_1_variable = aa_table_group %>% 
    group_by(specimen_id) %>% 
    filter(sum(n_aa == 1) == (loci_in_group - 1))
  
  # check if the above check found any groups with this scenario and if so then process 
  if(nrow(aa_table_group_only_1_variable) > 0){
    # first filter to the just the variable loci and calculate the within sample allele frequency (wsaf) 
    aa_table_group_only_1_variable_filt_variable = aa_table_group_only_1_variable %>% 
      filter(n_aa != 1) %>% 
      group_by(specimen_id, gene, gene_id, aa_position) %>% 
      mutate(total_read_count = sum(read_count)) %>%
      mutate(wsaf = read_count/total_read_count) %>% 
      group_by(specimen_id) %>% 
      mutate(within_sample_hap = row_number())
    
    # now grab the invariable calls to join with the variable loci, give it hap ID 
    # so it can be summarized later to create the full multi-locus call 
    aa_table_group_only_1_variable_filt_invariable = aa_table_group_only_1_variable %>% 
      filter(n_aa == 1) %>% 
      ungroup() %>% 
      select(specimen_id, gene_id, aa_position, aa, group_id) %>% 
      left_join(
        aa_table_group_only_1_variable_filt_variable %>% 
          ungroup() %>% 
          group_by(specimen_id) %>% 
          summarise(within_sample_hap = max(within_sample_hap)), 
        by = c("specimen_id")
      ) %>% 
      rowwise() %>% 
      mutate(within_sample_hap = list(1:within_sample_hap)) %>% 
      unnest(within_sample_hap)
    
    # combine with the variable 
    aa_table_group_only_1_variable_filt_combined = bind_rows(
      aa_table_group_only_1_variable_filt_variable, 
      aa_table_group_only_1_variable_filt_invariable
    )
    
    # now collapse the invariable and variable loci into the multi-locus calls 
    aa_table_group_only_1_variable_filt_combined_haps = aa_table_group_only_1_variable_filt_combined %>% 
      arrange(gene_id, aa_position, within_sample_hap) %>% 
      group_by(specimen_id, gene_id, within_sample_hap) %>% 
      summarise(positions = paste0(aa_position, collapse = "_"), 
                aas = paste0(aa, collapse = "_"), 
                wsaf = ifelse(all(is.na(wsaf)), NA, min(wsaf, na.rm = T))) %>% 
      unite(per_gene_variant, gene_id, positions, aas, sep = ":") %>% 
      group_by(specimen_id, within_sample_hap) %>% 
      summarise(variant = paste0(per_gene_variant, collapse = ";"), 
                wsaf = min(wsaf, na.rm = T))
  } else { 
    # if no specimens with only 1 variable call, create an empty tibble to join with the below 
    aa_table_group_only_1_variable_filt_combined_haps = tibble()
  }

  # filter to specimens that didn't have only 1 loci variable 
  aa_table_group_filt = aa_table_group %>% 
    filter(specimen_id %!in% aa_table_group_only_1_variable$specimen_id)
  
  # now we will determine multi-locus for the specimens by filtering to a specific within sample frequency 
  # and presume that the multi-locus haplotype exist if all sites have 1 call above this cut off 
  # cut off should be reasonably high to avoid false haplotype creation for example if setting to 0.51,
  # would have a high chance of creating false haplotypes
  # this filter will also capture samples that are completely monoclonal 
  aa_table_group_filt_dominant = aa_table_group_filt %>% 
    group_by(specimen_id, gene, gene_id, aa_position) %>% 
    mutate(total_read_count = sum(read_count)) %>%
    mutate(wsaf = read_count/total_read_count) %>% 
    filter(wsaf >= args$wsaf_cut_off) %>% 
    group_by(specimen_id) %>% 
    mutate(loci_called = n_distinct(paste0(gene_id, "-", aa_position))) %>% 
    group_by(specimen_id, gene, gene_id, aa_position) %>% 
    mutate(n_aa = n_distinct(aa)) %>% 
    filter(all(n_aa == 1), loci_called == loci_groups_split[[loci_group]]$loci_in_group[1])
  
  # collapse the calls to create the multi-locus call 
  aa_table_group_filt_dominant_collapse = aa_table_group_filt_dominant %>% 
    arrange(gene_id, aa_position) %>% 
    group_by(specimen_id, gene_id) %>% 
    summarise(positions = paste0(aa_position, collapse = "_"), 
              aas = paste0(aa, collapse = "_"), 
              wsaf = min(wsaf)) %>% 
    unite(per_gene_variant, gene_id, positions, aas, sep = ":") %>% 
    group_by(specimen_id) %>% 
    summarise(variant = paste0(per_gene_variant, collapse = ";"), 
              wsaf = min(wsaf)) %>% 
    group_by(specimen_id) %>% 
    mutate(within_sample_hap = row_number())
  
  # combine the two approaches to get a final list of multi-locus calls 
  aa_table_group_filt_final = bind_rows(aa_table_group_filt_dominant_collapse, 
                                        aa_table_group_only_1_variable_filt_combined_haps)
  
  # collect for all groups into one master table 
  all_aa_table_group_filt_final = bind_rows(
    all_aa_table_group_filt_final, 
    aa_table_group_filt_final %>% 
      mutate(group_id = loci_group)
  )
  
  # now calculate the sample prevalence and allele frequencies
  aa_table_group_filt_final_prev_freq <- switch(
    args$method,
    wsaf_prop = calculate_multilocus_af_prev_wsaf_prop(aa_table_group_filt_final),
    presence_absence = calculate_multilocus_af_prev_presence_absence(aa_table_group_filt_final)
  ) %>% 
    mutate(group_id = loci_group) 
    
  # combine into final table for all groups 
  all_aa_table_group_filt_final_prev_freq  = bind_rows(
    all_aa_table_group_filt_final_prev_freq, 
    aa_table_group_filt_final_prev_freq
  )
}

# write out the results 
write_tsv(all_aa_table_group_filt_final_prev_freq, args$output_path)

if(!is.null(args$recalc_single_locus_output_path)){
  # if an output path is set for export, also re-calculate single loci frequencies from the multi-locus calls 
  # these re-calculated frequencies can serve as a sanity check against the single locus freqs/prevs calculated directly
  # from the data, for example if some allele frequencies aren't present in the re-calcuated calls but are in the direct calculations 
  # then the multi-locus processing failed to capture any multi-locus haplotypes with that allele and is not capture the full data 
  slaf_from_ml = 
  slaf_from_ml <- switch(
    args$method,
    wsaf_prop = generate_single_locus_prev_freq_from_multilocus_groups_wsaf_prop(all_aa_table_group_filt_final),
    presence_absence = generate_single_locus_prev_freq_from_multilocus_groups_presence_absence(all_aa_table_group_filt_final)
  )
  write_tsv(slaf_from_ml, args$recalc_single_locus_output_path)
}


