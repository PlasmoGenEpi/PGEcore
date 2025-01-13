
# Example of estimate_coi_naive module

#------------------------------------------------------

library(optparse)

# Use default parameters from the options list
arg <- parse_args(OptionParser(option_list = opts))

# Manually add or edit some arguments
arg$input_path <- "././data/example2_allele_table.tsv"
arg$method <- "quantile_method"

# Estimate COI
df_coi <- run_estimate_coi_naive(input_path = arg$input_path,
                                 method = arg$method,
                                 integer_threshold = arg$integer_threshold,
                                 quantile_threshold = arg$quantile_threshold)


# Write to file
write_coi_naive(df_coi = df_coi,
                output_path = "././data/naive_coi_table.tsv")
