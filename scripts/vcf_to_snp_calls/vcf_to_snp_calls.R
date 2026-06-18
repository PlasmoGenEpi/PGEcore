#!/usr/bin/env Rscript
# vcf_to_snp_calls.R
# -----------------------------------------------------------------------------
# Reverse of snp_calls_to_vcf.R: read a VCF and emit per-sample SNP calls in a
# format as close as possible to pileup_specific_snps.R's snp_calls.tsv.gz. Read
# counts are sourced from FORMAT/AD, one row per (specimen, SNP, observed allele)
# where AD >= --min_reads.
#
# Only *simple* SNP alleles are emitted (single-base REF and ALT in {A,C,G,T});
# MNPs, indels, and symbolic alleles (e.g. *, <NON_REF>) are skipped. With
# --biallelic, only sites with exactly one simple ALT (REF + 1 ALT) are kept.
#
# Columns that have no source in a VCF (seq, he, etc.) are dropped: a VCF does
# not carry microhaplotype context. `strand` is always "+" because REF/ALT are
# written on the forward (genome) strand. `target_name` is emitted only if the
# variant's INFO carries a TARGET= key, otherwise the column is omitted.
# -----------------------------------------------------------------------------

packagesToLoad = c("tibble", "dplyr", "stringr", "readr", "optparse")
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

#' TRUE if x is a single canonical DNA base (a simple-SNP allele)
is_simple_base <- function(x){
  grepl("^[ACGTacgt]$", x)
}

#' Extract the TARGET= value from a VCF INFO string, or NA if absent
info_target <- function(info){
  m <- regmatches(info, regexec("(?:^|;)TARGET=([^;]+)", info))[[1]]
  if(length(m) >= 2) m[2] else NA_character_
}

# Parse arguments ------------------------------------------------------
opts <- list(
  make_option(
    "--vcf",
    help = "Input VCF (.vcf or .vcf.gz). Must carry FORMAT/AD (allelic depths)."
  ),
  make_option(
    "--snp_calls_output",
    help = "Output SNP-calls TSV path; gzip-compressed if it ends in .gz"
  ),
  make_option(
    "--biallelic",
    action = "store_true", default = FALSE,
    help = "Keep only biallelic sites (REF + exactly one simple ALT); multi-allelic sites are skipped"
  ),
  make_option(
    "--min_reads",
    type = "integer", default = 1L,
    help = "Minimum AD reads for an allele to be emitted as an observed call [default %default]"
  ),
  make_option(
    "--overwrite",
    action = "store_true", default = FALSE,
    help = "Overwrite --snp_calls_output if it already exists"
  ),
  make_option(
    "--verbose",
    action = "store_true", default = FALSE,
    help = "Print a summary message when finished (otherwise run silently)"
  )
)

run_vcf_to_snp_calls <- function(){
  arg <- parse_args(OptionParser(option_list = opts))
  checkOptparseRequiredArgsThrow(arg, c("vcf", "snp_calls_output"))

  if(file.exists(arg$snp_calls_output) && !arg$overwrite){
    stop(paste0(arg$snp_calls_output, " already exists, use --overwrite to overwrite"))
  }

  # Read the VCF (gzip-transparent) and split off the header.
  con <- if(grepl("\\.gz$", arg$vcf)) gzfile(arg$vcf, "rt") else file(arg$vcf, "rt")
  lines <- readLines(con)
  close(con)

  header_idx <- which(grepl("^#CHROM", lines))
  if(length(header_idx) == 0){
    stop("no #CHROM header line found in ", arg$vcf)
  }
  header <- strsplit(sub("^#", "", lines[header_idx[1]]), "\t")[[1]]
  if(length(header) < 10){
    stop("VCF has no sample columns; nothing to emit (need FORMAT + at least one sample)")
  }
  samples <- as.character(header[10:length(header)])

  data_lines <- lines[(header_idx[1] + 1):length(lines)]
  data_lines <- data_lines[nzchar(data_lines)]

  n_simple_sites <- 0L
  n_skipped_complex <- 0L
  n_skipped_multiallelic <- 0L
  chunks <- vector("list", length(data_lines))

  for(i in seq_along(data_lines)){
    f <- strsplit(data_lines[i], "\t")[[1]]
    s_chrom <- f[1]; s_pos <- as.integer(f[2]); s_id <- f[3]
    s_ref <- f[4]; s_alt <- f[5]; s_info <- f[8]; s_format <- f[9]
    gts <- f[10:length(f)]

    # require a single-base REF; otherwise this is an MNP/indel/symbolic site
    if(!is_simple_base(s_ref)){ n_skipped_complex <- n_skipped_complex + 1L; next }

    alts <- strsplit(s_alt, ",")[[1]]
    alt_simple <- vapply(alts, is_simple_base, logical(1))
    n_simple_alt <- sum(alt_simple)
    if(n_simple_alt == 0){ n_skipped_complex <- n_skipped_complex + 1L; next }
    if(arg$biallelic && n_simple_alt != 1){
      n_skipped_multiallelic <- n_skipped_multiallelic + 1L; next
    }

    # allele bases indexed 0=REF, 1..=ALT; keep REF + simple ALTs only
    allele_bases <- c(s_ref, alts)
    keep_allele  <- c(TRUE, alt_simple)               # always keep REF
    site_biallelic <- (n_simple_alt == 1)
    target_val <- info_target(s_info)

    # AD position within the per-sample FORMAT fields
    fmt <- strsplit(s_format, ":")[[1]]
    ad_idx <- match("AD", fmt)
    if(is.na(ad_idx)){ stop("FORMAT has no AD field at ", s_chrom, ":", s_pos) }

    # per-sample x per-allele depth matrix
    ad_mat <- matrix(NA_integer_, nrow = length(samples), ncol = length(allele_bases))
    for(j in seq_along(samples)){
      sub <- strsplit(gts[j], ":")[[1]]
      ad <- if(length(sub) >= ad_idx) sub[ad_idx] else "."
      if(is.na(ad) || ad == "." || ad == ""){ next }
      vals <- suppressWarnings(as.integer(strsplit(ad, ",")[[1]]))
      n <- min(length(vals), ncol(ad_mat))
      if(n > 0){ ad_mat[j, seq_len(n)] <- vals[seq_len(n)] }
    }

    # 0-based start / 1-based end, mirroring pileup snp_name when ID is absent
    pos0 <- s_pos - 1L
    snp_name <- if(is.na(s_id) || s_id == ".") sprintf("%s-%d-%d", s_chrom, pos0, s_pos) else s_id

    # emit a row per (sample, kept allele) where AD >= min_reads
    rows <- vector("list", sum(keep_allele))
    r <- 0L
    for(a in which(keep_allele)){
      depth <- ad_mat[, a]
      hit <- which(!is.na(depth) & depth >= arg$min_reads)
      if(length(hit) == 0){ next }
      r <- r + 1L
      rows[[r]] <- tibble(
        specimen_name = samples[hit],
        target_name   = target_val,
        chrom         = s_chrom,
        pos           = pos0,
        snp_name      = snp_name,
        strand        = "+",
        ref_base      = toupper(s_ref),
        seq_base      = toupper(allele_bases[a]),
        reads         = depth[hit],
        is_biallelic  = site_biallelic
      )
    }
    if(r > 0){
      chunks[[i]] <- bind_rows(rows[seq_len(r)])
      n_simple_sites <- n_simple_sites + 1L
    }
  }

  snp_calls <- bind_rows(chunks) |>
    arrange(specimen_name, chrom, pos, seq_base)

  # Drop target_name entirely unless at least one variant carried INFO TARGET=
  if("target_name" %in% colnames(snp_calls) && all(is.na(snp_calls$target_name))){
    snp_calls <- snp_calls |> select(-target_name)
  }

  readr::write_tsv(snp_calls, arg$snp_calls_output)

  if(arg$verbose){
    message(sprintf(
      "Wrote %d SNP calls (%d simple sites, %d samples) to %s (skipped %d non-SNP%s)",
      nrow(snp_calls), n_simple_sites, length(samples), arg$snp_calls_output,
      n_skipped_complex,
      if(arg$biallelic) sprintf(", %d multi-allelic", n_skipped_multiallelic) else ""
    ))
  }
  return(TRUE)
}

run_status <- run_vcf_to_snp_calls()
