library(dplyr)
library(tidyr)
library(readr)
library(optparse)
library(stringr)
library(tibble)

alfy_create_SNP_slice_input <- function(input_path, output_ref, output_alt) {
    # modified from a chatGPT conversion of a python script I wrote
    # Read the input file
    data <- read_tsv(input_path, show_col_types = FALSE)
    # Create gene_position identifier
    data <- data %>% 
      mutate(gene_pos = paste(gene_id, aa_position, sep = "_"))
    # Identify valid sites with fewer than 3 unique alleles
    good_alleles <- data %>% 
      group_by(gene_pos) %>% 
      summarize(unique_alleles = n_distinct(aa)) %>% 
      filter(unique_alleles < 3) %>% 
      pull(gene_pos)

    # Filter data to include only good alleles
    data <- data %>% 
      filter(gene_pos %in% good_alleles)

    # Compute ref and alt counts
    ref_alt_counts <- data %>% 
      mutate(ref_count = ifelse(ref_aa == aa, read_count, 0),
             alt_count = ifelse(ref_aa != aa, read_count, 0)) %>%
      group_by(specimen_id, gene_pos) %>% 
      summarize(ref_count = sum(ref_count),
                alt_count = sum(alt_count), .groups = "drop")
    

    # Pivot to wide format
    ref_counts <- ref_alt_counts %>% 
      select(specimen_id, gene_pos, ref_count) %>% 
      pivot_wider(names_from = gene_pos, values_from = ref_count, values_fill = 0)

    alt_counts <- ref_alt_counts %>% 
      select(specimen_id, gene_pos, alt_count) %>% 
      pivot_wider(names_from = gene_pos, values_from = alt_count, values_fill = 0)

    # Write to files
    write.table(ref_counts, output_ref, sep = "\t", row.names = FALSE, quote = FALSE)
    write.table(alt_counts, output_alt, sep = "\t", row.names = FALSE, quote = FALSE)
}

alfy_create_ref_alt <-function(input_path, amino_acids) {
    # Read the input file
    data <- read_tsv(input_path, show_col_types = FALSE)
    ref_alt <- data %>%
      mutate(gene_pos = paste(gene_id, aa_position, sep = "_")) %>%
      distinct(gene_pos, ref_aa, aa) %>%
      filter(ref_aa != aa)
    write_tsv(ref_alt, amino_acids)
}


run_SNPslice_neg_binomial <- function(input_data_path1, input_data_path2, output_dir, gap, nmcmc = 10000, alpha = 2, rep = 1, script_path = "my_snpslicemain.R") {
    #borrowed from Kathryn
    model <- 3
    # tmp_file_path <- file.path(output_dir, "tmp_formatted_output.txt")
    # # If model is not 0 then need read counts, so need to alter this function to take a second input path
    # create_SNPslice_input(input_data_path, tmp_file_path)
    # Run snp slice
    command <- paste(
        "Rscript ", script_path,
        " model=", model, " nmcmc=", nmcmc, " alpha=", alpha, " gap=", gap,
        " input=", input_data_path1, " input2=", input_data_path2, " output_dir=", output_dir,
        sep = ""
    )
    print(command)
    runtime <- system.time({
        system(command)
    })
    # file.remove(tmp_file_path)
    # Infer the output
    D_filename <- paste("neg_D_nmcmc", nmcmc, "_gap", gap, "_rep", rep, ".txt", sep = "")
    A_filename <- paste("neg_A_nmcmc", nmcmc, "_gap", gap, "_rep", rep, ".txt", sep = "")
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
    print('********************************')
    size <- nrow(collapsed_table)
    sum_column <- rep(total_sum, size)
    print('size is')
    print(size)
    new_collapsed <- collapsed_table %>% add_column(sum_column)
    print('new collapsed is')
    print(new_collapsed)
    return(new_collapsed)
}


# Function to get alternates
get_alternates <- function(assignments, alt_counts) {
  alt_dict <- read_tsv(assignments, col_names = c("gene_pos", "ref_aa", "alt_aa")) %>%
    select(gene_pos, ref_aa, alt_aa)
  
  loci <- read_tsv(alt_counts, n_max = 1) %>%
    colnames() %>%
    .[-1]
  
  loci_df <- tibble(locus = loci) %>%
    left_join(alt_dict, by = c("locus" = "gene_pos"))
  
  return(loci_df)
}

# Function to convert genotype to STAVE format
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

# Function to lookup genotype
lookup_genotype <- function(genotype, genotype_mappings) {
  tibble(
    char = unlist(str_split(genotype, "")),
    locus = genotype_mappings$locus
  ) %>%
    mutate(
      alt = ifelse(char == "0", genotype_mappings$ref_aa, genotype_mappings$alt_aa)
    ) %>%
    select(locus, alt)
}


# Process frequencies and generate output
genotypes_from_freqs <- function(freq_estimates, output_file, genotype_mappings) {
  freqs <- read_tsv(freq_estimates)
  
  header_dict <- names(freqs)
  
  output <- freqs %>%
    rowwise() %>%
    mutate(
      found_genotype = list(lookup_genotype(plsf_table.sequence, genotype_mappings)),
      stave_string = genotype_to_stave(found_genotype)
    ) %>%
    ungroup() %>%
    select(stave_string, plsf_table.SNPslice_frequency, plsf_table.sum_column) %>%
    rename(
      frequency = plsf_table.SNPslice_frequency,
      total_genotype_count = plsf_table.sum_column
    )
  
  write_tsv(output, output_file)
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
    "--output", 
    help = "final stave output path"
  ),



  make_option(
    "--gap", 
    help = "how many iterations of no improvement before stopping"
  )
)

arg <- parse_args(OptionParser(option_list = opts))
print('args are')
print(arg)
assignments <- "alt_alleles.tsv"

alfy_create_ref_alt(arg$aa_calls, assignments)
output_ref='ref_counts_R_v2.tsv'
output_alt='alt_counts_R_v2.tsv'
alfy_create_SNP_slice_input(arg$aa_calls, output_ref, output_alt)
output_list=run_SNPslice_neg_binomial(output_alt, output_ref, getwd(), gap=arg$gap)
freq=output_list[1]
freq_estimates <- 'snp-slice_freqs.tsv'

write.table(freq, freq_estimates, sep = "\t", row.names = FALSE, quote = FALSE)

# File paths (update these accordingly)
output_file <- arg$output

# Read genotype mappings
genotype_mappings <- get_alternates(assignments, output_alt)
# Execute the function
genotypes_from_freqs(freq_estimates, output_file, genotype_mappings)