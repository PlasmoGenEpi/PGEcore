library(dplyr)
library(tidyr)
library(readr)
library(optparse)
library(stringr)
library(tibble)
library(roxygen2)

#' Create SNP-Slice Input Data Frames
#'
#' Processes amino acid calls to generate reference and alternative read count matrices
#' in the format required by SNP-Slice.
#'
#' @param data A data frame with at least the columns:
#'   \code{specimen_id}, \code{gene_id}, \code{aa_position}, \code{aa},
#'   \code{ref_aa}, and \code{read_count}.
#'
#' @return A list containing:
#'   \item{ref_counts}{Wide data frame of reference allele counts per specimen.}
#'   \item{alt_counts}{Wide data frame of alternative allele counts per specimen.}
#' @export
create_SNP_slice_input_df <- function(data) {
  data <- data %>% mutate(gene_pos = paste(gene_id, aa_position, sep = "_"))

  good_alleles <- data %>%
    group_by(gene_pos) %>%
    summarize(unique_alleles = n_distinct(aa), .groups = "drop") %>%
    filter(unique_alleles < 3) %>%
    pull(gene_pos)

  data <- data %>% filter(gene_pos %in% good_alleles)

  ref_alt_counts <- data %>%
    mutate(
      ref_count = ifelse(ref_aa == aa, read_count, 0),
      alt_count = ifelse(ref_aa != aa, read_count, 0)
    ) %>%
    group_by(specimen_id, gene_pos) %>%
    summarize(across(c(ref_count, alt_count), sum), .groups = "drop")

  list(
    ref_counts = ref_alt_counts %>%
      select(specimen_id, gene_pos, ref_count) %>%
      pivot_wider(names_from = gene_pos, values_from = ref_count, values_fill = 0),
    alt_counts = ref_alt_counts %>%
      select(specimen_id, gene_pos, alt_count) %>%
      pivot_wider(names_from = gene_pos, values_from = alt_count, values_fill = 0)
  )
}


#' Create Reference/Alternate Allele Mapping
#'
#' Extracts a unique set of loci with their reference and alternate amino acids.
#'
#' @param data Data frame with \code{gene_id}, \code{aa_position}, \code{ref_aa}, and \code{aa}.
#'
#' @return A data frame with columns \code{gene_pos}, \code{ref_aa}, and \code{alt_aa}.
#' @export
create_ref_alt_df <- function(data) {
  data %>%
    mutate(gene_pos = paste(gene_id, aa_position, sep = "_")) %>%
    distinct(gene_pos, ref_aa, aa) %>%
    filter(ref_aa != aa) %>%
    rename(alt_aa = aa)
}

#' Run SNP-Slice
#'
#' Executes the SNP-Slice R script with specified parameters and parses results.
#'
#' @param gap Integer; number of no-improvement iterations before stopping.
#' @param script_path Path to the adapted SNP-Slice main R script.
#' @param ref Path to the reference allele counts input file.
#' @param alt Path to the alternative allele counts input file.
#' @param snp_slice_dir Path to the cloned SNP-Slice GitHub repository.
#' @param model Model number to use in SNP-Slice.
#' @param nmcmc Number of MCMC iterations.
#' @param alpha Alpha parameter for SNP-Slice.
#' @param rep Replicate number.
#'
#' @return A list containing:
#'   \item{plsf_table}{Data frame of collapsed SNP-Slice frequencies.}
#'   \item{runtime}{Execution time.}
#' @export
run_SNPslice <- function(gap, script_path, ref, alt, snp_slice_dir, model = 3, nmcmc = 10000, alpha = 2, rep = 1) {
  output_dir <- "output" # this output folder is hardcoded in snpslicemain.R so needs to be hardcoded as the same value here
  dir.create(file.path(output_dir), showWarnings = FALSE)
  command <- paste(
    "Rscript ", script_path,
    " model=", model, " nmcmc=", nmcmc, " alpha=", alpha, " gap=", gap,
    " input=", ref, " input2=", alt, " output_dir=", output_dir, " snp_slice_code=", snp_slice_dir,
    sep = ""
  )
  runtime <- system.time({
    system(command)
  })
  D_filename <- paste(output_dir, "/neg_D_nmcmc", nmcmc, "_gap", gap, "_rep", rep, ".txt", sep = "")
  A_filename <- paste(output_dir, "/neg_A_nmcmc", nmcmc, "_gap", gap, "_rep", rep, ".txt", sep = "")
  haplotype_dict_path <- D_filename
  host_strain_association <- A_filename
  plsf_table <- infer_SNPslice_freqs(haplotype_dict_path, host_strain_association)
  return(
    list(
      plsf_table = plsf_table,
      runtime = runtime
    )
  )
}

#' Infer SNP-Slice Frequencies
#'
#' Reads SNP-Slice haplotype and host-strain association outputs, and calculates
#' haplotype frequencies.
#'
#' @param haplotype_dict Path to the SNP-Slice haplotype dictionary output file.
#' @param host_strain_association Path to the SNP-Slice host-strain association output file.
#'
#' @return A data frame with columns \code{sequence}, \code{SNPslice_frequency},
#'   and \code{sum_column}.
#' @export
infer_SNPslice_freqs <- function(haplotype_dict, host_strain_association) {
  # Read allele association output
  D <- read.table(haplotype_dict, header = FALSE, sep = "\t")
  A <- read.table(host_strain_association, header = FALSE, sep = "\t")
  # Calculate the frequency for each of the haplotypes
  column_sums <- colSums(A)
  total_sum <- sum(A)
  column_sum_ratios <- column_sums / total_sum
  row_strings <- apply(D, 1, paste, collapse = "")
  combined_table <- cbind(row_strings, column_sum_ratios)
  colnames(combined_table) <- c("sequence", "SNPslice_frequency")
  combined_table <- as.data.frame(combined_table, stringsAsFactors = FALSE)
  combined_table$SNPslice_frequency <- as.numeric(combined_table$SNPslice_frequency)
  # Collapse the table by summing SNPslice_frequency for each sequence
  collapsed_table <- aggregate(SNPslice_frequency ~ sequence, combined_table, sum)
  size <- nrow(collapsed_table)
  sum_column <- rep(total_sum, size)
  new_collapsed <- collapsed_table %>% add_column(sum_column)
  return(new_collapsed)
}


#' Get Alternate Alleles Data Frame
#'
#' Aligns loci names from the SNP-Slice alt counts table with allele assignments.
#'
#' @param assignments_df Data frame from \code{create_ref_alt_df()}.
#' @param alt_counts_df Wide alternative counts table (first column \code{specimen_id}).
#'
#' @return A data frame mapping loci to reference and alternate alleles.
#' @export
get_alternates_df <- function(assignments_df, alt_counts_df) {
  loci <- colnames(alt_counts_df)[-1] # skip specimen_id
  tibble(locus = loci) %>%
    left_join(assignments_df, by = c("locus" = "gene_pos"))
}




#' Convert Genotype to Stave String
#'
#' Orders genotype loci and produces a concise stave-format string.
#'
#' @param found_genotype Data frame with columns \code{locus} and \code{alt}.
#'
#' @return A character string encoding the genotype in stave format.
#' @export
genotype_to_stave <- function(found_genotype) {
  found_genotype %>%
    mutate(
      gene = str_extract(locus, "^.*(?=_[0-9]+$)"),
      pos = as.integer(str_extract(locus, "[0-9]+$"))
    ) %>%
    arrange(gene, pos) %>%
    group_by(gene) %>%
    summarize(
      stave_string = paste(
        str_c(pos, collapse = "_"),
        str_c(alt, collapse = ""),
        sep = ":"
      )
    ) %>%
    summarize(stave_string = str_c(stave_string, collapse = ";")) %>%
    pull(stave_string)
}

#' Lookup Genotype Alleles
#'
#' Translates a binary genotype string (0 = reference, 1 = alternate)
#' into a table of loci and amino acids.
#'
#' @param genotype Character string of 0/1 alleles.
#' @param genotype_mappings Data frame mapping loci to \code{ref_aa} and \code{alt_aa}.
#'
#' @return A data frame with columns \code{locus} and \code{alt}.
#' @export
lookup_genotype <- function(genotype, genotype_mappings) {
  tibble::tibble(
    char = unlist(stringr::str_split(genotype, "")),
    locus = genotype_mappings$locus
  ) %>%
    dplyr::mutate(
      alt = ifelse(char == "0", genotype_mappings$ref_aa, genotype_mappings$alt_aa)
    ) %>%
    dplyr::select(locus, alt)
}


#' Generate Stave Genotypes from SNP-Slice Frequencies
#'
#' Maps SNP-Slice output sequences to their corresponding stave string format.
#'
#' @param freq_table Data frame from \code{infer_SNPslice_freqs()}.
#' @param genotype_mappings Data frame from \code{get_alternates_df()}.
#'
#' @return A data frame with columns \code{stave_string}, \code{frequency},
#'   and \code{total_genotype_count}.
#' @export
genotypes_from_freqs_df <- function(freq_table, genotype_mappings) {
  freq_table %>%
    rowwise() %>%
    mutate(
      found_genotype = list(lookup_genotype(sequence, genotype_mappings)),
      stave_string = genotype_to_stave(found_genotype)
    ) %>%
    ungroup() %>%
    select(stave_string, SNPslice_frequency, sum_column) %>%
    rename(
      frequency = SNPslice_frequency,
      total_genotype_count = sum_column
    )
}


#' Process and Subset Loci Groups
#'
#' Splits amino acid call data into loci groups, runs SNP-Slice on each group,
#' and combines results into a single output file.
#'
#' @param aa_calls_path Path to TSV of amino acid calls.
#' @param loci_group_table_path Path to TSV of loci group assignments.
#'
#' @return None; writes combined results to the file specified in \code{arg$output}.
#' @export
subset_groups <- function(aa_calls_path, loci_group_table_path) {
  dir.create(file.path("inputdata"), showWarnings = FALSE)

  aa_calls <- read_tsv(aa_calls_path, show_col_types = FALSE)
  loci_groups <- read_tsv(loci_group_table_path, show_col_types = FALSE)

  results_list <- list()

  for (group in unique(loci_groups$group_id)) {
    message("Processing group: ", group)

    # Filter for this group
    loci_sub <- loci_groups %>% filter(group_id == group)
    filtered_aa_calls <- aa_calls %>%
      semi_join(loci_sub, by = c("gene_id", "aa_position"))

    # Create assignments in memory
    assignments_df <- create_ref_alt_df(filtered_aa_calls)

    # Create SNP-Slice ref/alt tables in memory
    snp_input <- create_SNP_slice_input_df(filtered_aa_calls)

    # Write only what SNP-Slice needs
    ref_path <- file.path("inputdata", paste0(group, "_ref.txt"))
    alt_path <- file.path("inputdata", paste0(group, "_alt.txt"))
    write.table(snp_input$ref_counts, ref_path, sep = "\t", row.names = FALSE, quote = FALSE)
    write.table(snp_input$alt_counts, alt_path, sep = "\t", row.names = FALSE, quote = FALSE)

    # Run SNP-Slice
    output_list <- run_SNPslice(
      gap = arg$gap,
      script_path = "scripts/snp-slice_wrapper/adapted_snpslicemain.R",
      ref = ref_path,
      alt = alt_path,
      snp_slice_dir = arg$snp_slice_dir,
      model = arg$model,
      rep = 1
    )

    # Map alternates (in memory)
    genotype_mappings <- get_alternates_df(assignments_df, snp_input$alt_counts)

    # Convert to stave (in memory)
    group_results <- genotypes_from_freqs_df(output_list$plsf_table, genotype_mappings) %>%
      mutate(group_id = group)

    results_list[[group]] <- group_results
  }

  # Combine & write once at the end
  combined <- bind_rows(results_list)
  write_tsv(combined, arg$output)
  message("All groups processed and combined output saved to ", arg$output)
}


opts <- list(
  make_option(
    "--aa_calls",
    help = str_c(
      "TSV containing amino acid calls, with the columns: specimen_id, ",
      "target_id, gene_id, aa_position, ref_codon, ref_aa, codon, aa"
    )
  ),
  make_option(
    "--loci_group_table",
    type = "character",
    help = str_c(
      "TSV containing loci groups, with the columns: group_id, gene_id,
      aa_position."
    )
  ),
  make_option(
    "--snp_slice_dir",
    help = "path to cloned snp-slice github repo"
  ),
  make_option(
    "--output",
    help = "final stave output path"
  ),
  make_option(
    "--model",
    help = "which model number to use, see snp-slice github for documentation"
  ),
  make_option(
    "--gap",
    help = "how many iterations of no improvement before stopping"
  )
)

working_directory <- getwd()
arg <- parse_args(OptionParser(option_list = opts))

# Run snp-slice per group
subset_groups(arg$aa_calls, arg$loci_group_table)
