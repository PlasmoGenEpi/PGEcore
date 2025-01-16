#!/usr/bin/env Rscript


packagesToLoad = c("tibble", "dplyr", "stringr", "readr", "optparse", "pwalign")

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

#' get the aligned pos for the real position of the sequence, adjusting for gaps
#'
#' @param aligned_seq the aligned seq with gaps
#' @param pos the "real" position
#'
#' @returns the align position in the aligned_seq that corresponds to the real position
getAlnPosPerRealPos <- function(aligned_seq, pos){
  gapless_pos <- which(strsplit(as.character(aligned_seq), NULL)[[1]] != "-")
  gapless_pos[pos]
}


#' For overlap type alignment get the aligned bases of the pattern including the end gaps
#'
#' @param pw_overlapAlign the overlap alignment
#'
#' @returns the aligned pattern with end gaps
getAlignedPatternFromOverlapAlign <-function(pw_overlapAlign){
  overlapAlign_alignedPattern =  alignedPattern(pw_overlapAlign)
  
  overlapAlign_alignedPattern_bases = as.character(overlapAlign_alignedPattern)
  
  preceding_gap_size = pw_overlapAlign@subject@range@start - 1 
  trailing_gap_size = nchar(pw_overlapAlign@subject@unaligned) - (pw_overlapAlign@subject@range@start + pw_overlapAlign@subject@range@width) + 1
  
  preceding_other_gap_size = pw_overlapAlign@pattern@range@start - 1 
  trailing_other_gap_size = nchar(pw_overlapAlign@pattern@unaligned) - (pw_overlapAlign@pattern@range@start + pw_overlapAlign@pattern@range@width) + 1
  
  if(preceding_other_gap_size > 0){
    overlapAlign_alignedPattern_bases = paste0(substr(as.character(pw_overlapAlign@pattern@unaligned), 1, preceding_other_gap_size), overlapAlign_alignedPattern_bases)
  }
  if(trailing_other_gap_size > 0){
    overlapAlign_alignedPattern_bases = paste0(overlapAlign_alignedPattern_bases, substr(as.character(pw_overlapAlign@pattern@unaligned), nchar(pw_overlapAlign@pattern@unaligned) - trailing_other_gap_size + 1, nchar(pw_overlapAlign@pattern@unaligned)))
  }
  paste0(paste0(rep("-", preceding_gap_size),collapse = ""), overlapAlign_alignedPattern_bases, paste0(rep("-", trailing_gap_size), collapse = ""))
}


#' For overlap type alignment get the aligned bases of the subject including the end gaps
#'
#' @param pw_overlapAlign the overlap alignment 
#'
#' @returns the aligned subject with end gaps 
getAlignedSubjectFromOverlapAlign <-function(pw_overlapAlign){
  overlapAlign_alignedSubject =  alignedSubject(pw_overlapAlign)
  
  overlapAlign_alignedSubject_bases = as.character(overlapAlign_alignedSubject)
  
  preceding_gap_size = pw_overlapAlign@pattern@range@start - 1 
  trailing_gap_size = nchar(pw_overlapAlign@pattern@unaligned) - (pw_overlapAlign@pattern@range@start + pw_overlapAlign@pattern@range@width) + 1
  
  preceding_other_gap_size = pw_overlapAlign@subject@range@start - 1 
  trailing_other_gap_size = nchar(pw_overlapAlign@subject@unaligned) - (pw_overlapAlign@subject@range@start + pw_overlapAlign@subject@range@width) + 1
  
  if(preceding_other_gap_size > 0){
    overlapAlign_alignedSubject_bases = paste0(substr(as.character(pw_overlapAlign@subject@unaligned), 1, preceding_other_gap_size), overlapAlign_alignedSubject_bases)
  }
  if(trailing_other_gap_size > 0){
    overlapAlign_alignedSubject_bases = paste0(overlapAlign_alignedSubject_bases, substr(as.character(pw_overlapAlign@subject@unaligned), nchar(pw_overlapAlign@subject@unaligned) - trailing_other_gap_size + 1, nchar(pw_overlapAlign@subject@unaligned)))
  }
  paste0(paste0(rep("-", preceding_gap_size),collapse = ""), overlapAlign_alignedSubject_bases, paste0(rep("-", trailing_gap_size), collapse = ""))
}


# Parse arguments ------------------------------------------------------
opts <- list(
  make_option(
    "--allele_table", 
    help = str_c(
      "TSV containing the columns: specimen_id, target_id, read_count, seq"
    )
  ),
  make_option(
    "--ref_bed", 
    help = str_c(
      "a bed file containing the reference location of the ref_seq, no column names but the first 6 columns should be chrom, start, end, target_id, length, strand, ref_seq"
    )
  ), 
  make_option(
    "--snps_of_interest", 
    help = str_c(
      "a bed file containing the snps of interest to calculate pileup, is genomic location of the snp position so all loci should be of size 1, should have columns, #chrom, start, end, name, length, strand"
    )
  ), 
  make_option(
    "--output_directory", 
    help = str_c(
      "the output directory to write results to"
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
    "--overwrite_dir", 
    help = str_c(
      "overwrite the output directory if it already exists"
    ), 
    action="store_true", 
    default=FALSE,
    type = "logical"
  ), 
  make_option(
    "--collapse_calls_by_summing", 
    help = str_c(
      "collapse amino acid calls by summing across targets, by default the target with the highest read count will be used as the final call"
    ), 
    action="store_true", 
    default=FALSE,
    type = "logical"
  )
)

required_arguments = c("allele_table", "ref_bed", "snps_of_interest", "output_directory")

# parse arguments
arg <- parse_args(OptionParser(option_list = opts))

## check for required arguments
checkOptparseRequiredArgsThrow(arg, required_arguments)

# create output directory and decided to overwrite or not if it exists 
if(dir.exists(arg$output_directory) & arg$overwrite_dir){
  unlink(arg$output_directory, recursive = T)
} else if(dir.exists(arg$output_directory)){
  stop(paste0(arg$output_directory, " already exist, use --overwrite_dir to overwrite"))
}
dir.create(arg$output_directory)

# process if only processing specific targets and specimens

select_target_ids = c()
if("select_target_ids" %in% names(arg)){
  if(file.exists(arg$select_target_ids)){
    select_target_ids = readr::read_tsv(arg$select_target_ids, col_names = F)$X1
  } else { 
    select_target_ids = unlist(strsplit(arg$select_target_ids, split = ","))
  }
}


select_specimen_ids = c()
if("select_specimen_ids" %in% names(arg)){
  if(file.exists(arg$select_specimen_ids)){
    select_specimen_ids = readr::read_tsv(arg$select_specimen_ids, col_names = F)$X1
  } else { 
    select_specimen_ids = unlist(strsplit(arg$select_specimen_ids, split = ","))
  }
}


# read in input data and gather warnings about columns and data formatting 
warnings = c()

# read in the panel reference information 
ref_bed = readr::read_tsv(arg$ref_bed, col_names = T)
warnings = c(warnings, genWarningsMissCols(ref_bed, c("#chrom", "start", "end", "target_id", "length", "strand", "ref_seq"), arg$ref_bed))

if(length(select_target_ids) > 0 ){
  missing_sel_tars = setdiff(select_target_ids, ref_bed$target_id)
  if(length(missing_sel_tars) > 0 ){
    warnings = c(warnings, paste0("supplied --select_target_ids but the following targets are missing from ", arg$ref_bed, "\n", 
                                  paste0(missing_sel_tars, collapse = ",")))
  }
  ref_bed = ref_bed |> 
    filter(target_id %in% select_target_ids)
}

# read in the loci of interest 
snps_of_interest = readr::read_tsv(arg$snps_of_interest, col_names = T)
warnings = c(warnings, genWarningsMissCols(snps_of_interest, c("#chrom", "start", "end", "name", "length", "strand"), arg$snps_of_interest))

## check to make sure loci are of length 3 only 
for(row in 1:nrow(snps_of_interest)){
  if((snps_of_interest$end[row] - snps_of_interest$start[row]) != 1){
    warnings  = c(warnings, paste0("loci of interest must be of length 1, locus: ", snps_of_interest$name[row], " is length: ",  snps_of_interest$length[row]))
  }
}

# read in allele table for the microhaplotype data 
allele_table = readr::read_tsv(arg$allele_table)
warnings = c(warnings, genWarningsMissCols(allele_table, c("specimen_id","target_id","read_count","seq"), arg$allele_table))

if(length(select_specimen_ids) > 0 ){
  missing_sel_specs = setdiff(select_target_ids, unique(allele_table$specimen_id))
  if(length(missing_sel_specs) > 0 ){
    warnings = c(warnings, paste0("supplied --select_specimen_ids but the following specimen_ids are missing from ", arg$allele_table, "\n", 
                                  paste0(missing_sel_specs, collapse = ",")))
  }
  allele_table = allele_table |> 
    filter(specimen_id %in% select_specimen_ids)
}

if(length(select_target_ids) > 0 ){
  missing_sel_tars = setdiff(select_target_ids, unique(allele_table$target_id))
  if(length(missing_sel_tars) > 0 ){
    warnings = c(warnings, paste0("supplied --select_target_ids but the following target_ids are missing from ", arg$allele_table, "\n", 
                                  paste0(missing_sel_tars, collapse = ",")))
  }
  allele_table = allele_table |> 
    filter(target_id %in% select_target_ids)
}


# check to see if values between dataets are similar 
ref_allele_decomp = set_decompose(ref_bed$target_id, unique(allele_table$target_id))

if(length(ref_allele_decomp$only_in_vectorB) > 0){
  warnings = c(warnings, paste0("the following loci were missing from the reference location file ", arg$ref_bed, " but are in ", arg$allele_table,
                                "\n", paste0(ref_allele_decomp$only_in_vectorB, collapse = ",")
  ) 
  )
}

if(length(warnings) > 0){
  stop(paste0("\n", paste0(warnings, collapse = "\n")) )
}

# create a table of unique 
allele_table_unique_haps = allele_table |> 
  select(target_id, seq) |> 
  unique() |> 
  arrange(target_id)


# find which loci and targets intersect
snps_of_interest_by_target_loci = list()

ref_bed$intersected_snps_of_interest = ""
snps_of_interest$covered_by_target = ""

for(row in 1:nrow(ref_bed)){
  for(loci_row in 1:nrow(snps_of_interest)){
    if(ref_bed$`#chrom`[row] == snps_of_interest$`#chrom`[loci_row] & 
       snps_of_interest$start[loci_row] >= ref_bed$start[row] & 
       snps_of_interest$end[loci_row] <= ref_bed$end[row]){
      # add which loci are covered by the target to the target table
      if ("" != ref_bed$intersected_snps_of_interest[row]){
        ref_bed$intersected_snps_of_interest[row] = paste0(ref_bed$intersected_snps_of_interest[row], ",")
      }
      ref_bed$intersected_snps_of_interest[row] = paste0(ref_bed$intersected_snps_of_interest[row], loci_row)
      
      # add which targets cover this loci
      if ("" != snps_of_interest$covered_by_target[loci_row]){
        snps_of_interest$covered_by_target[loci_row] = paste0(snps_of_interest$covered_by_target[loci_row], ",")
      }
      snps_of_interest$covered_by_target[loci_row] = paste0(snps_of_interest$covered_by_target[loci_row], ref_bed$target_id[row])
    }
  }
}
snps_of_interest = snps_of_interest |> 
  mutate(covered_by_target = ifelse("" == covered_by_target, "uncovered", covered_by_target))

# create a map of target location to key into with target IDs 
ref_bed_by_loci = list()
for(row in 1:nrow(ref_bed)){
  ref_bed_by_loci[[ref_bed$target_id[row]]] = ref_bed[row,]
}


ref_bed_withInterest = ref_bed |> 
  filter("" != intersected_snps_of_interest)

microhaps_with_snps_of_interest = ref_bed_withInterest$target_id

# create substitution matrix 
mat <- pwalign::nucleotideSubstitutionMatrix(match = 2, mismatch = -2, baseOnly = TRUE)



all_snps_of_interest_for_target_for_microhap = tibble()

# iterate over unique sequences to translate
for(row in 1:nrow(allele_table_unique_haps)){
  # if the target is in the loci of interest table then determine the calls per sequence 
  if(allele_table_unique_haps$target_id[row] %in% microhaps_with_snps_of_interest){
    # end-to-end align sequences 
    overlapAlign <- pwalign::pairwiseAlignment(Biostrings::DNAString(allele_table_unique_haps$seq[row]), 
                                               Biostrings::DNAString(ref_bed_by_loci[[allele_table_unique_haps$target_id[row]]]$ref_seq[1]),
                                               substitutionMatrix = mat, gapOpening = 5, gapExtension = 1, 
                                               type="overlap") 
    # get the loci for this target 
    snps_of_interest_for_target = snps_of_interest[as.numeric(unlist(strsplit(ref_bed_by_loci[[allele_table_unique_haps$target_id[row]]]$intersected_snps_of_interest, ","))),] |> 
      mutate(rel_start = start - ref_bed_by_loci[[allele_table_unique_haps$target_id[row]]]$start[1], 
             rel_end = end - ref_bed_by_loci[[allele_table_unique_haps$target_id[row]]]$start[1])
    
    snps_of_interest_for_target_for_microhap = tibble()
    # get the relative position of the codon within the aligned sequence and translate 
    for(snps_of_interest_for_target_row in 1:nrow(snps_of_interest_for_target)){
      aln_pos = getAlnPosPerRealPos(getAlignedSubjectFromOverlapAlign(overlapAlign), snps_of_interest_for_target$rel_start[snps_of_interest_for_target_row] + 1)
      seq_base = Biostrings::DNAString(substr(getAlignedPatternFromOverlapAlign(overlapAlign), aln_pos, aln_pos + snps_of_interest_for_target$length[snps_of_interest_for_target_row] -1))
      ref_base = Biostrings::DNAString(substr(ref_bed_by_loci[[allele_table_unique_haps$target_id[row]]]$ref_seq[1], 
                                               snps_of_interest_for_target$rel_start[snps_of_interest_for_target_row] + 1, 
                                               snps_of_interest_for_target$rel_start[snps_of_interest_for_target_row] + 1 + snps_of_interest_for_target$length[snps_of_interest_for_target_row] - 1))
      # if the sequences are in the opposite direction from the amino acid call then reverse complement to get the correct translation 
      if(snps_of_interest_for_target$strand[snps_of_interest_for_target_row] != ref_bed_by_loci[[allele_table_unique_haps$target_id[row]]]$strand){
        seq_base = Biostrings::reverseComplement(seq_base)
        ref_base = Biostrings::reverseComplement(ref_base)
      }

      # create the table with the data of interest 
      snps_of_interest_for_target_for_microhap = 
        bind_rows(
          snps_of_interest_for_target_for_microhap,
          tibble(
            target_id = allele_table_unique_haps$target_id[row], 
            seq = allele_table_unique_haps$seq[row], 
            chrom = snps_of_interest_for_target$`#chrom`[snps_of_interest_for_target_row], 
            pos = snps_of_interest_for_target$start[snps_of_interest_for_target_row],
            snp_name = snps_of_interest_for_target$name[snps_of_interest_for_target_row],
            ref_base = as.character(ref_base), 
            seq_base = as.character(seq_base)
          )
        )
    }
    
    # join in with the final summary table 
    all_snps_of_interest_for_target_for_microhap = bind_rows(
      all_snps_of_interest_for_target_for_microhap, 
      snps_of_interest_for_target_for_microhap
    )
  }
}


# take the calls per unique sequence and join them to the original allele calls for all samples 
allele_table_out = allele_table |> 
  filter(target_id %in% microhaps_with_snps_of_interest) |> 
  left_join(all_snps_of_interest_for_target_for_microhap, relationship = "many-to-many", by = c("target_id", "seq"))

# get sample coverage info 
coveredBySamplesCount = allele_table_out |> 
  group_by(chrom, pos, snp_name, ref_base) |>
  summarise(n_samples = n_distinct(specimen_id)) |> 
  mutate(total_samples = n_distinct(allele_table$specimen_id))

snps_of_interest_out = snps_of_interest |> 
  left_join(coveredBySamplesCount |>
              dplyr::rename(name = snp_name), 
            by = c("name"))


allele_table_out_uncallable = allele_table_out |> 
  filter(seq_base == "-")

allele_table_out_filt = allele_table_out |> 
  filter(seq_base != "-")

arg$collapse_calls_by_summing =F
# collapse amino acid calls 
if(arg$collapse_calls_by_summing){
  allele_table_out_collapsed = allele_table_out_filt |> 
    group_by(specimen_id, chrom, pos, snp_name, ref_base, seq_base) |> 
    summarise(read_count = sum(read_count), 
              target_id = paste0(target_id, collapse = ","))
} else { 
  allele_table_out_winnerTarget = allele_table_out_filt |> 
    group_by(specimen_id, chrom, pos, snp_name, ref_base, target_id) |> 
    summarise(read_count = sum(read_count)) |> 
    arrange(desc(read_count)) |> 
    mutate(read_count_rank = row_number(), 
           covered_by_target_ids = paste0(target_id, collapse = ",")) |> 
    filter(read_count_rank == 1) |> 
    ungroup() |> 
    select(-read_count_rank) |> 
    dplyr::rename(best_target_id = target_id)
  
  allele_table_out_collapsed = allele_table_out_filt |> 
    left_join(allele_table_out_winnerTarget |> 
                ungroup() |> 
                select(-read_count), 
              by = c("specimen_id", "chrom", "pos", "snp_name", "ref_base")) |> 
    filter(target_id == best_target_id) |> 
    select(-seq)
  
  allele_table_out_collapsed = allele_table_out_collapsed |> 
    group_by(specimen_id, target_id, chrom, pos, snp_name, ref_base, seq_base, best_target_id, covered_by_target_ids) |> 
    summarise(read_count = sum(read_count))
}


# adding biallelic info 

allele_table_out = allele_table_out %>% 
  group_by(chrom, pos, snp_name, ref_base) %>% 
  mutate(isBiallelic = n_distinct(seq_base) <=2)

allele_table_out_collapsed = allele_table_out_collapsed %>% 
  group_by(chrom, pos, snp_name, ref_base) %>% 
  mutate(isBiallelic = n_distinct(seq_base) <=2)

# writing output results 
write_tsv(allele_table_out, file.path(arg$output_directory, "snp_calls.tsv.gz"))

# writing output results collapsing calls if calls occur on multiple targets  
write_tsv(allele_table_out_collapsed, file.path(arg$output_directory,"collapsed_snp_calls.tsv.gz"))

# writing covered by info 
write_tsv(snps_of_interest_out, file.path(arg$output_directory,"loci_covered_by_target_samples_info.tsv"))

if(nrow(allele_table_out_uncallable) > 0){
  write_tsv(allele_table_out_uncallable, file.path(arg$output_directory,"allele_table_out_uncallable.tsv"))
}

