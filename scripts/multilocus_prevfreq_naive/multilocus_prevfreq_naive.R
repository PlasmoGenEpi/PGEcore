#!/usr/bin/env Rscript 

# quietly load packages 
packagesToLoad = c("dplyr", "optparse", "purrr", "readr", "stringr", "tidyr", "validate")

loaded = suppressMessages(lapply(packagesToLoad, require, character.only = TRUE))

options(readr.show_col_types = FALSE)
options(dplyr.summarise.inform = FALSE)

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
    is.character(gene_id), 
    is.integer(aa_position), 
    is.integer(read_count), 
    is.character(aa), 
    ! is.na(specimen_id), 
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

  # tidy up columns
  # df_aa <- df_aa |>
  #   select(specimen_id, gene_id, aa_position, read_count, aa) #|>
    # rename(gene = gene_id,
    #        pos = aa_position)
  
  # get numer of amino acids at each locus
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


#' Generate the single locus allele frequencies and prevalence from a table of multilocus calls 
#'
#' @param multilocus_calls A table of 3 columns, group_id, specimen_id, and variant (in variantstring format)
#'
#' @returns summarized per group, single locus freq/prev recalculated from the multilocus variant calls 
#' @export
#'
#' @examples
generate_sinlge_locus_prev_freq_from_multilocus_groups <-function(multilocus_calls){
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
      "Input multilocus_calls for function generate_sinlge_locus_prev_freq_from_multilocus_groups failed one or more validation checks: ", 
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

# RUN MODULE -----------------------------------------------------------

# parse arguments
args <- parse_args(OptionParser(option_list = opts))
# Arguments used for development
if (interactive()) {
  args = list()
  args$aa_table <- "../../data/example_amino_acid_calls.tsv"
  args$loci_groups_input <- "../../data/example_loci_groups.tsv"
  args$output_path <- "../../mlafp.tsv"
}

# args = list()
# args$aa_table <- "/Users/nicholashathaway/Dropbox (Personal)/ownCloud/documents/plasmodium/falciparum/ucsf/uganda_jessica/IMMERSE_Post_Analysis/R1_R2_R3_R4_ASV_data_thomas_subsetted_filtered_for_analysis_1B_2_only_2025_09_26_unfiltered/work/de/568f81412cd776e4fe96ff1e19e6e1/Atiak-2023.collapsed_amino_acid_calls.tsv.gz"
# args$loci_groups_input <- "/Users/nicholashathaway/Dropbox (Personal)/ownCloud/documents/plasmodium/falciparum/ucsf/uganda_jessica/IMMERSE_Post_Analysis/R1_R2_R3_R4_ASV_data_thomas_subsetted_filtered_for_analysis_1B_2_only_2025_09_26_unfiltered/crt_dhfr_dhps_loci_groups.tsv"
# args$output_path <- "~/Documents/sourceCodes/plasmodiumdrugres/bin/PGEcore/scripts/multilocus_prevfreq_naive/testing_mlafp_more_groups2.tsv"
# args$recalc_single_locus_output_path <- "~/Documents/sourceCodes/plasmodiumdrugres/bin/PGEcore/scripts/multilocus_prevfreq_naive/testing_slaf_from_mlafp_more_groups2.tsv"
# args$wsaf_cut_off = 0.70

required_arguments = c("aa_table", "loci_groups_input", "output_path")

checkOptparseRequiredArgsThrow(args, required_arguments)

# Read in data
aa_table <- create_aa_table_input(args$aa_table)


loci_groups <- create_loci_group_input(args$loci_groups_input) %>% 
  group_by(group_id) %>% 
  mutate(loci_in_group = n_distinct(paste0(gene_id, "-", aa_position)))

loci_groups_split = split(loci_groups, loci_groups$group_id)

all_aa_table_group_filt_final_prev_freq = tibble()

all_aa_table_group_filt_final = tibble()


for(loci_group in names(loci_groups_split)){
  # filter to loci in group and filter to only samples that have calls for all loci in the group 
  aa_table_group = aa_table %>% 
    inner_join(loci_groups_split[[loci_group]], by = c("gene_id", "aa_position")) %>% 
    group_by(specimen_id) %>% 
    mutate(loci_called = n_distinct(paste0(gene_id, "-", aa_position))) %>% 
    filter(loci_called == loci_in_group) %>% 
    ungroup()
  
  aa_table_group_only_1_variable = aa_table_group %>% 
    group_by(specimen_id) %>% 
    filter(sum(n_aa == 1) == (loci_in_group - 1))
  if(nrow(aa_table_group_only_1_variable) > 0){
    aa_table_group_only_1_variable_filt_variable = aa_table_group_only_1_variable %>% 
      filter(n_aa != 1) %>% 
      group_by(specimen_id, gene, gene_id, aa_position) %>% 
      mutate(total_read_count = sum(read_count)) %>%
      mutate(wsaf = read_count/total_read_count) %>% 
      group_by(specimen_id) %>% 
      mutate(within_sample_hap = row_number())
    
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
    
    aa_table_group_only_1_variable_filt_combined = bind_rows(
      aa_table_group_only_1_variable_filt_variable, 
      aa_table_group_only_1_variable_filt_invariable
    )
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
    aa_table_group_only_1_variable_filt_combined_haps = tibble()
  }

  aa_table_group_filt = aa_table_group %>% 
    filter(specimen_id %!in% aa_table_group_only_1_variable$specimen_id)
  
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
  
  aa_table_group_filt_dominant_collapse = aa_table_group_filt_dominant %>% 
    arrange(gene_id, aa_position) %>% 
    group_by(specimen_id, gene_id) %>% 
    summarise(positions = paste0(aa_position, collapse = "_"), 
              aas = paste0(aa, collapse = "_")) %>% 
    unite(per_gene_variant, gene_id, positions, aas, sep = ":") %>% 
    group_by(specimen_id) %>% 
    summarise(variant = paste0(per_gene_variant, collapse = ";")) %>% 
    group_by(specimen_id) %>% 
    mutate(within_sample_hap = row_number()) %>% 
    mutate(wsaf = 1)
  
  aa_table_group_filt_final = bind_rows(aa_table_group_filt_dominant_collapse, 
                                        aa_table_group_only_1_variable_filt_combined_haps)
  
  all_aa_table_group_filt_final = bind_rows(
    all_aa_table_group_filt_final, 
    aa_table_group_filt_final %>% 
      mutate(group_id = loci_group)
  )
  
  aa_table_group_filt_final_prev_freq = aa_table_group_filt_final %>% 
    ungroup() %>% 
    mutate(sample_total = n_distinct(specimen_id), 
           allele_total = n()) %>% 
    group_by(variant, sample_total, allele_total) %>% 
    summarise(allele_count = n(), 
           sample_count = n_distinct(specimen_id)) %>% 
    mutate(prev = sample_count/sample_total, 
           freq = allele_count/allele_total) %>% 
    ungroup() %>% 
    mutate(group_id = loci_group) %>% 
    select(group_id, variant, allele_count, sample_count, allele_total, sample_total, freq, prev)
    
  all_aa_table_group_filt_final_prev_freq  = bind_rows(
    all_aa_table_group_filt_final_prev_freq, 
    aa_table_group_filt_final_prev_freq
  )
}

write_tsv(all_aa_table_group_filt_final_prev_freq, args$output_path)

if(!is.null(args$recalc_single_locus_output_path)){
  slaf_from_ml = generate_sinlge_locus_prev_freq_from_multilocus_groups(all_aa_table_group_filt_final)
  write_tsv(slaf_from_ml, args$recalc_single_locus_output_path)
}


