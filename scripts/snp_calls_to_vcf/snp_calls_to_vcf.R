#!/usr/bin/env Rscript
# snp_calls_to_vcf.R
# -----------------------------------------------------------------------------
# Build a VCF from pileup SNP calls produced by pileup_specific_snps.R. The VCF
# carries per-sample allelic depths (FORMAT/AD) so downstream tools (e.g.
# moimix::getFws) can compute within-host statistics.
#
# Works on EITHER the raw `snp_calls.tsv.gz` or the `collapsed_snp_calls.tsv.gz`
# output: it only needs the columns common to both (specimen_name, chrom, pos,
# snp_name, ref_base, seq_base, reads) and sums reads per (specimen, SNP, allele)
# across overlapping targets -- a real sum for the raw file, a no-op for the
# already-collapsed file.
#
# Multi-allelic sites are emitted with comma-separated ALT and multi-value AD
# (use --biallelic to keep only sites with a single ALT). GT is derived from AD
# at the requested --ploidy. The reference genome FASTA supplies accurate
# ##contig lengths.
# -----------------------------------------------------------------------------

packagesToLoad = c("tibble", "dplyr", "stringr", "readr", "optparse", "Biostrings")
loaded = suppressMessages(lapply(packagesToLoad, require, character.only = TRUE))
options(readr.show_col_types = FALSE)
options(dplyr.summarise.inform = FALSE)

#' Stop if any required args are missing
checkOptparseRequiredArgsThrow <- function(arg, required_args){
  missing <- setdiff(required_args, names(arg))
  if(length(missing) > 0){
    stop(paste0("missing the following arguments: ", paste0("--", missing, collapse = ", ")))
  }
}

#' Columns in `columns` not present in `tib`
returnMissingColumns <- function(tib, columns){ setdiff(columns, colnames(tib)) }

#' Reverse-complement a vector of DNA strings (handles single bases)
revcomp <- function(x){
  comp <- chartr("ACGTacgtNn", "TGCAtgcaNn", x)
  vapply(strsplit(comp, NULL), function(ch) paste(rev(ch), collapse = ""), character(1))
}

#' Derive a fixed-ploidy GT string from per-allele depths.
#'
#' Alleles with depth >= min_reads are "present". They fill the `ploidy` slots
#' ordered by depth: more present than ploidy -> keep the top `ploidy`; fewer ->
#' pad with the most-supported present allele; none present -> all-missing (./.).
#' Returned indices are 0-based (0 = REF) and sorted, joined by "/".
derive_gt <- function(depths, ploidy, min_reads){
  present <- which(depths >= min_reads)                 # 1-based positions
  if(length(present) == 0){
    return(paste(rep(".", ploidy), collapse = "/"))
  }
  present <- present[order(depths[present], decreasing = TRUE)]
  if(length(present) >= ploidy){
    chosen <- present[1:ploidy]
  } else {
    chosen <- c(present, rep(present[1], ploidy - length(present)))
  }
  paste(sort(chosen - 1L), collapse = "/")
}

# Parse arguments ------------------------------------------------------
opts <- list(
  make_option(
    "--snp_calls",
    help = str_c(
      "TSV of SNP calls from pileup_specific_snps.R (raw snp_calls.tsv.gz or ",
      "collapsed_snp_calls.tsv.gz). Required columns: specimen_name, chrom, pos, ",
      "snp_name, strand, ref_base, seq_base, reads"
    )
  ),
  make_option(
    "--genome",
    help = "Reference genome FASTA; used to write accurate ##contig=<ID=,length=> headers"
  ),
  make_option(
    "--vcf_output",
    help = "Output VCF path; gzip-compressed if it ends in .gz"
  ),
  make_option(
    "--biallelic",
    action = "store_true", default = FALSE,
    help = "Keep only biallelic sites (REF + exactly one ALT); multi-allelic sites are skipped"
  ),
  make_option(
    "--ploidy",
    type = "integer", default = 2L,
    help = "Ploidy used to render the GT field [default %default]"
  ),
  make_option(
    "--gt_min_reads",
    type = "integer", default = 1L,
    help = "Minimum reads for an allele to count as present when deriving GT [default %default]"
  ),
  make_option(
    "--overwrite",
    action = "store_true", default = FALSE,
    help = "Overwrite --vcf_output if it already exists"
  ),
  make_option(
    "--verbose",
    action = "store_true", default = FALSE,
    help = "Print a summary message when finished (otherwise run silently)"
  )
)

run_snp_calls_to_vcf <- function(){
  arg <- parse_args(OptionParser(option_list = opts))
  checkOptparseRequiredArgsThrow(arg, c("snp_calls", "genome", "vcf_output"))

  if(file.exists(arg$vcf_output) && !arg$overwrite){
    stop(paste0(arg$vcf_output, " already exists, use --overwrite to overwrite"))
  }

  calls <- readr::read_tsv(arg$snp_calls, col_types = readr::cols(specimen_name = readr::col_character()))
  needed <- c("specimen_name", "chrom", "pos", "snp_name", "strand", "ref_base", "seq_base", "reads")
  miss <- returnMissingColumns(calls, needed)
  if(length(miss) > 0){
    stop(paste0("snp_calls is missing required columns: ", paste0(miss, collapse = ", ")))
  }

  # contig names + lengths from the genome (names truncated at first whitespace)
  genome <- Biostrings::readDNAStringSet(arg$genome)
  genome_names <- sub("\\s.*$", "", names(genome))
  genome_lengths <- setNames(Biostrings::width(genome), genome_names)

  # drop uncallable ("-") / missing-allele rows, then sum reads per
  # (specimen, SNP, allele) across overlapping targets to get allele depths.
  # (a real sum for the raw file; a no-op for the already-collapsed file)
  ad <- calls |>
    filter(!is.na(seq_base), seq_base != "-", !is.na(reads)) |>
    group_by(specimen_name, chrom, pos, snp_name, strand, ref_base, seq_base) |>
    summarise(reads = sum(reads), .groups = "drop")

  samples <- sort(unique(ad$specimen_name))
  snps <- ad |>
    distinct(chrom, pos, snp_name, strand, ref_base) |>
    arrange(chrom, pos)

  contigs <- sort(unique(ad$chrom))
  missing_contigs <- setdiff(contigs, genome_names)
  if(length(missing_contigs) > 0){
    stop(paste0("contigs in snp_calls not found in --genome: ", paste0(missing_contigs, collapse = ", ")))
  }

  # VCF header
  header <- c(
    "##fileformat=VCFv4.2",
    "##source=PGEcore snp_calls_to_vcf.R",
    paste0("##contig=<ID=", contigs, ",length=", genome_lengths[contigs], ">"),
    "##INFO=<ID=NS,Number=1,Type=Integer,Description=\"Number of samples with data (depth > 0)\">",
    "##INFO=<ID=DP,Number=1,Type=Integer,Description=\"Total read depth across samples\">",
    "##INFO=<ID=AF,Number=A,Type=Float,Description=\"Alt allele frequency from pooled read depths\">",
    "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype derived from allelic depths\">",
    "##FORMAT=<ID=AD,Number=R,Type=Integer,Description=\"Allelic depths for the ref and alt alleles in the order listed\">",
    "##FORMAT=<ID=DP,Number=1,Type=Integer,Description=\"Total read depth\">",
    paste0("#", paste(c("CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT", samples), collapse = "\t"))
  )

  lines <- character(nrow(snps))
  n_lines <- 0L
  n_monomorphic <- 0L
  n_multiallelic_skipped <- 0L
  missing_field <- paste0(paste(rep(".", arg$ploidy), collapse = "/"), ":.:.")

  for(i in 1:nrow(snps)){
    s_chrom <- snps$chrom[i]; s_pos <- snps$pos[i]
    s_name <- snps$snp_name[i]; s_ref <- snps$ref_base[i]; s_strand <- snps$strand[i]

    sub <- ad[ad$chrom == s_chrom & ad$pos == s_pos &
              ad$snp_name == s_name & ad$ref_base == s_ref, ]
    alts <- sort(setdiff(unique(sub$seq_base), s_ref))
    if(length(alts) == 0){ n_monomorphic <- n_monomorphic + 1L; next }   # skip monomorphic
    if(arg$biallelic && length(alts) > 1){
      n_multiallelic_skipped <- n_multiallelic_skipped + 1L; next
    }
    # bases from the pileup are in the SNP's strand orientation; emit forward
    # (genome) strand so REF/ALT match the reference. Depths are looked up with
    # the strand-oriented bases (`alleles`), so AD order still matches REF,ALT.
    alleles <- c(s_ref, alts)
    if(s_strand == "-"){
      ref_out <- revcomp(s_ref); alts_out <- revcomp(alts)
    } else {
      ref_out <- s_ref; alts_out <- alts
    }

    # depth vector per sample (NULL if the sample has no reads at this site)
    sub_by_spec <- split(sub, sub$specimen_name)
    depth_list <- lapply(samples, function(sm){
      srows <- sub_by_spec[[sm]]
      if(is.null(srows)){ return(NULL) }
      as.integer(vapply(alleles, function(a) sum(srows$reads[srows$seq_base == a]), numeric(1)))
    })
    names(depth_list) <- samples

    # site-level INFO from pooled depths
    present <- depth_list[!vapply(depth_list, is.null, logical(1))]
    allele_totals <- if(length(present) > 0) Reduce(`+`, present) else integer(length(alleles))
    site_dp <- sum(allele_totals)
    alt_totals <- allele_totals[-1]
    af <- if(site_dp > 0) alt_totals / site_dp else rep(0, length(alt_totals))
    info <- sprintf("NS=%d;DP=%d;AF=%s", length(present), site_dp,
                    paste(formatC(af, format = "g", digits = 6), collapse = ","))

    # per-sample FORMAT fields
    fields <- vapply(samples, function(sm){
      depths <- depth_list[[sm]]
      if(is.null(depths)){ return(missing_field) }
      paste0(derive_gt(depths, arg$ploidy, arg$gt_min_reads), ":",
             paste(depths, collapse = ","), ":", sum(depths))
    }, character(1))

    n_lines <- n_lines + 1L
    lines[n_lines] <- paste(c(s_chrom, as.integer(s_pos) + 1L, s_name, ref_out,
                              paste(alts_out, collapse = ","), ".", ".", info,
                              "GT:AD:DP", fields), collapse = "\t")
  }
  lines <- lines[seq_len(n_lines)]

  con <- if(grepl("\\.gz$", arg$vcf_output)) gzfile(arg$vcf_output, "w") else file(arg$vcf_output, "w")
  writeLines(c(header, lines), con)
  close(con)

  if(arg$verbose){
    message(sprintf(
      "Wrote %d variant sites for %d samples to %s (skipped %d monomorphic%s)",
      n_lines, length(samples), arg$vcf_output, n_monomorphic,
      if(arg$biallelic) sprintf(", %d multi-allelic", n_multiallelic_skipped) else ""
    ))
  }
  return(TRUE)
}

run_status <- run_snp_calls_to_vcf()
