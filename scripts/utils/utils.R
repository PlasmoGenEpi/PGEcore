suppressPackageStartupMessages(library(dplyr))

convert_to_stave <- function(df) {

  df %>%
    ungroup %>%
    mutate(variant = paste(gene_id, aa_position, aa, sep = ":")) %>%
    select(variant, prev)

}