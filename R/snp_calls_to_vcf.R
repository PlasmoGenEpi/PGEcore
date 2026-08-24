#' Reverse-complement a vector of DNA strings (handles single bases)
#'
#' @param x Character vector of DNA sequences.
#' @return Character vector of reverse complements.
#' @keywords internal
revcomp <- function(x) {
  comp <- chartr("ACGTacgtNn", "TGCAtgcaNn", x)
  vapply(
    strsplit(comp, NULL),
    function(ch) paste(rev(ch), collapse = ""),
    character(1)
  )
}

#' Derive a fixed-ploidy GT string from per-allele depths
#'
#' Alleles with depth `>= min_reads` are "present". They fill the `ploidy` slots
#' ordered by depth: more present than ploidy keeps the top `ploidy`; fewer pads
#' with the most-supported present allele; none present yields all-missing
#' (`./.`). Returned indices are 0-based (`0` = REF) and sorted, joined by `/`.
#'
#' @param depths Integer vector of depths in REF, ALT order.
#' @param ploidy Ploidy used to render GT.
#' @param min_reads Minimum reads for an allele to count as present.
#' @return A GT string such as `"0/1"` or `"./."`.
#' @keywords internal
derive_gt <- function(depths, ploidy, min_reads) {
  present <- which(depths >= min_reads)
  if (length(present) == 0) {
    return(paste(rep(".", ploidy), collapse = "/"))
  }
  present <- present[order(depths[present], decreasing = TRUE)]
  if (length(present) >= ploidy) {
    chosen <- present[seq_len(ploidy)]
  } else {
    chosen <- c(present, rep(present[1], ploidy - length(present)))
  }
  paste(sort(chosen - 1L), collapse = "/")
}

#' Build a VCF from pileup SNP calls
#'
#' Converts pileup SNP calls (raw or collapsed) into a VCF with per-sample
#' allelic depths (`FORMAT/AD`). Reads are summed per `(specimen, SNP, allele)`
#' across overlapping targets. `REF`/`ALT` are written on the forward (genome)
#' strand. Monomorphic sites are skipped. Requires **Biostrings** (Suggests) to
#' read `--genome` contig lengths.
#'
#' @param snp_calls Path to a SNP-calls TSV, or a data frame, with columns
#'   `specimen_name`, `chrom`, `pos`, `snp_name`, `strand`, `ref_base`,
#'   `seq_base`, and `reads`. `pos` is 0-based (emitted as 1-based VCF `POS`).
#' @param genome Path to a reference genome FASTA used for `##contig` lengths.
#' @param vcf_output Output VCF path; gzip-compressed if it ends in `.gz`.
#' @param biallelic If `TRUE`, keep only sites with a single ALT.
#' @param ploidy Ploidy used to render the GT field.
#' @param gt_min_reads Minimum reads for an allele to count as present in GT.
#' @param overwrite If `FALSE` (default), refuse to overwrite `vcf_output`.
#' @param verbose If `TRUE`, print a summary message when finished.
#'
#' @return Invisibly returns `TRUE`.
#'
#' @export
snp_calls_to_vcf <- function(snp_calls,
                             genome,
                             vcf_output,
                             biallelic = FALSE,
                             ploidy = 2L,
                             gt_min_reads = 1L,
                             overwrite = FALSE,
                             verbose = FALSE) {
  options(dplyr.summarise.inform = FALSE)
  stop_if_output_exists(vcf_output, overwrite)
  check_suggested_pkg(
    "Biostrings",
    "reading the reference genome FASTA via snp_calls_to_vcf()"
  )

  needed <- c(
    "specimen_name", "chrom", "pos", "snp_name", "strand",
    "ref_base", "seq_base", "reads"
  )
  if (is.character(snp_calls) && length(snp_calls) == 1L) {
    if (!file.exists(snp_calls)) {
      stop(snp_calls, " does not exist", call. = FALSE)
    }
    calls <- readr::read_tsv(
      snp_calls,
      col_types = readr::cols(specimen_name = readr::col_character()),
      show_col_types = FALSE
    )
  } else if (is.data.frame(snp_calls)) {
    calls <- tibble::as_tibble(snp_calls)
  } else {
    stop("`snp_calls` must be a file path or a data frame.", call. = FALSE)
  }
  validate_required_columns(calls, needed, "snp_calls")

  genome_set <- read_genome_dna_string_set(
    genome,
    "reading the reference genome FASTA via snp_calls_to_vcf()"
  )
  genome_names <- names(genome_set)
  genome_lengths <- stats::setNames(Biostrings::width(genome_set), genome_names)

  ad <- calls |>
    dplyr::filter(
      !is.na(.data$seq_base),
      .data$seq_base != "-",
      !is.na(.data$reads)
    ) |>
    dplyr::group_by(
      .data$specimen_name,
      .data$chrom,
      .data$pos,
      .data$snp_name,
      .data$strand,
      .data$ref_base,
      .data$seq_base
    ) |>
    dplyr::summarise(reads = sum(.data$reads), .groups = "drop")

  samples <- sort(unique(ad$specimen_name))
  snps <- ad |>
    dplyr::distinct(
      .data$chrom, .data$pos, .data$snp_name, .data$strand, .data$ref_base
    ) |>
    dplyr::arrange(.data$chrom, .data$pos)

  contigs <- sort(unique(ad$chrom))
  missing_contigs <- setdiff(contigs, genome_names)
  if (length(missing_contigs) > 0) {
    stop(
      "contigs in snp_calls not found in --genome: ",
      paste(missing_contigs, collapse = ", "),
      call. = FALSE
    )
  }

  header <- c(
    "##fileformat=VCFv4.2",
    "##source=PGEcore snp_calls_to_vcf.R",
    paste0("##contig=<ID=", contigs, ",length=", genome_lengths[contigs], ">"),
    paste0(
      "##INFO=<ID=NS,Number=1,Type=Integer,",
      "Description=\"Number of samples with data (depth > 0)\">"
    ),
    paste0(
      "##INFO=<ID=DP,Number=1,Type=Integer,",
      "Description=\"Total read depth across samples\">"
    ),
    paste0(
      "##INFO=<ID=AF,Number=A,Type=Float,",
      "Description=\"Alt allele frequency from pooled read depths\">"
    ),
    paste0(
      "##FORMAT=<ID=GT,Number=1,Type=String,",
      "Description=\"Genotype derived from allelic depths\">"
    ),
    paste0(
      "##FORMAT=<ID=AD,Number=R,Type=Integer,",
      "Description=\"Allelic depths for the ref and alt alleles in the order listed\">"
    ),
    paste0(
      "##FORMAT=<ID=DP,Number=1,Type=Integer,",
      "Description=\"Total read depth\">"
    ),
    paste0(
      "#",
      paste(
        c(
          "CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO",
          "FORMAT", samples
        ),
        collapse = "\t"
      )
    )
  )

  lines <- character(nrow(snps))
  n_lines <- 0L
  n_monomorphic <- 0L
  n_multiallelic_skipped <- 0L
  missing_field <- paste0(paste(rep(".", ploidy), collapse = "/"), ":.:.")

  for (i in seq_len(nrow(snps))) {
    s_chrom <- snps$chrom[i]
    s_pos <- snps$pos[i]
    s_name <- snps$snp_name[i]
    s_ref <- snps$ref_base[i]
    s_strand <- snps$strand[i]

    sub <- ad[
      ad$chrom == s_chrom & ad$pos == s_pos &
        ad$snp_name == s_name & ad$ref_base == s_ref,
    ]
    alts <- sort(setdiff(unique(sub$seq_base), s_ref))
    if (length(alts) == 0) {
      n_monomorphic <- n_monomorphic + 1L
      next
    }
    if (isTRUE(biallelic) && length(alts) > 1) {
      n_multiallelic_skipped <- n_multiallelic_skipped + 1L
      next
    }
    alleles <- c(s_ref, alts)
    if (identical(s_strand, "-")) {
      ref_out <- revcomp(s_ref)
      alts_out <- revcomp(alts)
    } else {
      ref_out <- s_ref
      alts_out <- alts
    }

    sub_by_spec <- split(sub, sub$specimen_name)
    depth_list <- lapply(samples, function(sm) {
      srows <- sub_by_spec[[sm]]
      if (is.null(srows)) {
        return(NULL)
      }
      as.integer(vapply(
        alleles,
        function(a) sum(srows$reads[srows$seq_base == a]),
        numeric(1)
      ))
    })
    names(depth_list) <- samples

    present <- depth_list[!vapply(depth_list, is.null, logical(1))]
    allele_totals <- if (length(present) > 0) {
      Reduce(`+`, present)
    } else {
      integer(length(alleles))
    }
    site_dp <- sum(allele_totals)
    alt_totals <- allele_totals[-1]
    af <- if (site_dp > 0) alt_totals / site_dp else rep(0, length(alt_totals))
    info <- sprintf(
      "NS=%d;DP=%d;AF=%s",
      length(present),
      site_dp,
      paste(formatC(af, format = "g", digits = 6), collapse = ",")
    )

    fields <- vapply(samples, function(sm) {
      depths <- depth_list[[sm]]
      if (is.null(depths)) {
        return(missing_field)
      }
      paste0(
        derive_gt(depths, ploidy, gt_min_reads), ":",
        paste(depths, collapse = ","), ":", sum(depths)
      )
    }, character(1))

    n_lines <- n_lines + 1L
    lines[n_lines] <- paste(
      c(
        s_chrom, as.integer(s_pos) + 1L, s_name, ref_out,
        paste(alts_out, collapse = ","), ".", ".", info,
        "GT:AD:DP", fields
      ),
      collapse = "\t"
    )
  }
  lines <- lines[seq_len(n_lines)]

  con <- if (grepl("\\.gz$", vcf_output)) {
    gzfile(vcf_output, "w")
  } else {
    file(vcf_output, "w")
  }
  writeLines(c(header, lines), con)
  close(con)

  if (isTRUE(verbose)) {
    message(sprintf(
      "Wrote %d variant sites for %d samples to %s (skipped %d monomorphic%s)",
      n_lines, length(samples), vcf_output, n_monomorphic,
      if (isTRUE(biallelic)) {
        sprintf(", %d multi-allelic", n_multiallelic_skipped)
      } else {
        ""
      }
    ))
  }
  invisible(TRUE)
}
