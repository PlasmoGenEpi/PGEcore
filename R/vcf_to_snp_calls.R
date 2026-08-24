#' TRUE if `x` is a single canonical DNA base (a simple-SNP allele)
#'
#' @param x Character string.
#' @return Logical.
#' @keywords internal
is_simple_base <- function(x) {
  grepl("^[ACGTacgt]$", x)
}

#' Extract the TARGET= value from a VCF INFO string, or NA if absent
#'
#' @param info VCF INFO field.
#' @return Character scalar or `NA_character_`.
#' @keywords internal
info_target <- function(info) {
  m <- regmatches(info, regexec("(?:^|;)TARGET=([^;]+)", info))[[1]]
  if (length(m) >= 2) m[2] else NA_character_
}

#' Convert a VCF with FORMAT/AD into pileup-style SNP calls
#'
#' Reverse of [snp_calls_to_vcf()]: emits one row per `(specimen, SNP, observed
#' allele)` where allelic depth is at least `min_reads`. Only simple SNP alleles
#' (single-base REF and ALT in `{A,C,G,T}`) are kept. `pos` is 0-based. `strand`
#' is always `+`. `target_name` is included only when INFO carries `TARGET=`.
#'
#' @param vcf Path to a VCF (`.vcf` or `.vcf.gz`) with `FORMAT/AD`.
#' @param snp_calls_output Output TSV path; gzip-compressed if it ends in `.gz`.
#'   If `NULL`, results are returned without writing.
#' @param biallelic If `TRUE`, keep only sites with exactly one simple ALT.
#' @param min_reads Minimum AD reads for an allele to be emitted.
#' @param overwrite If `FALSE` (default), refuse to overwrite `snp_calls_output`.
#' @param verbose If `TRUE`, print a summary message when finished.
#'
#' @return A tibble of SNP calls. When `snp_calls_output` is set, the table is
#'   also written to that path.
#'
#' @export
vcf_to_snp_calls <- function(vcf,
                             snp_calls_output = NULL,
                             biallelic = FALSE,
                             min_reads = 1L,
                             overwrite = FALSE,
                             verbose = FALSE) {
  options(dplyr.summarise.inform = FALSE)
  stop_if_output_exists(snp_calls_output, overwrite)

  if (!is.character(vcf) || length(vcf) != 1L) {
    stop("`vcf` must be a path to a VCF file.", call. = FALSE)
  }
  if (!file.exists(vcf)) {
    stop(vcf, " does not exist", call. = FALSE)
  }

  con <- if (grepl("\\.gz$", vcf)) gzfile(vcf, "rt") else file(vcf, "rt")
  lines <- readLines(con)
  close(con)

  header_idx <- which(grepl("^#CHROM", lines))
  if (length(header_idx) == 0) {
    stop("no #CHROM header line found in ", vcf, call. = FALSE)
  }
  header <- strsplit(sub("^#", "", lines[header_idx[1]]), "\t")[[1]]
  if (length(header) < 10) {
    stop(
      "VCF has no sample columns; nothing to emit (need FORMAT + at least one sample)",
      call. = FALSE
    )
  }
  samples <- as.character(header[10:length(header)])

  if (header_idx[1] < length(lines)) {
    data_lines <- lines[(header_idx[1] + 1):length(lines)]
    data_lines <- data_lines[nzchar(data_lines)]
  } else {
    data_lines <- character(0)
  }

  n_simple_sites <- 0L
  n_skipped_complex <- 0L
  n_skipped_multiallelic <- 0L
  chunks <- vector("list", length(data_lines))

  for (i in seq_along(data_lines)) {
    f <- strsplit(data_lines[i], "\t")[[1]]
    s_chrom <- f[1]
    s_pos <- as.integer(f[2])
    s_id <- f[3]
    s_ref <- f[4]
    s_alt <- f[5]
    s_info <- f[8]
    s_format <- f[9]
    gts <- f[10:length(f)]

    if (!is_simple_base(s_ref)) {
      n_skipped_complex <- n_skipped_complex + 1L
      next
    }

    alts <- strsplit(s_alt, ",")[[1]]
    alt_simple <- vapply(alts, is_simple_base, logical(1))
    n_simple_alt <- sum(alt_simple)
    if (n_simple_alt == 0) {
      n_skipped_complex <- n_skipped_complex + 1L
      next
    }
    if (isTRUE(biallelic) && n_simple_alt != 1) {
      n_skipped_multiallelic <- n_skipped_multiallelic + 1L
      next
    }

    allele_bases <- c(s_ref, alts)
    keep_allele <- c(TRUE, alt_simple)
    site_biallelic <- (n_simple_alt == 1)
    target_val <- info_target(s_info)

    fmt <- strsplit(s_format, ":")[[1]]
    ad_idx <- match("AD", fmt)
    if (is.na(ad_idx)) {
      stop("FORMAT has no AD field at ", s_chrom, ":", s_pos, call. = FALSE)
    }

    ad_mat <- matrix(NA_integer_, nrow = length(samples), ncol = length(allele_bases))
    for (j in seq_along(samples)) {
      sub <- strsplit(gts[j], ":")[[1]]
      ad <- if (length(sub) >= ad_idx) sub[ad_idx] else "."
      if (is.na(ad) || ad == "." || ad == "") {
        next
      }
      vals <- suppressWarnings(as.integer(strsplit(ad, ",")[[1]]))
      n <- min(length(vals), ncol(ad_mat))
      if (n > 0) {
        ad_mat[j, seq_len(n)] <- vals[seq_len(n)]
      }
    }

    pos0 <- s_pos - 1L
    snp_name <- if (is.na(s_id) || s_id == ".") {
      sprintf("%s-%d-%d", s_chrom, pos0, s_pos)
    } else {
      s_id
    }

    rows <- vector("list", sum(keep_allele))
    r <- 0L
    for (a in which(keep_allele)) {
      depth <- ad_mat[, a]
      hit <- which(!is.na(depth) & depth >= min_reads)
      if (length(hit) == 0) {
        next
      }
      r <- r + 1L
      rows[[r]] <- tibble::tibble(
        specimen_name = samples[hit],
        target_name = target_val,
        chrom = s_chrom,
        pos = pos0,
        snp_name = snp_name,
        strand = "+",
        ref_base = toupper(s_ref),
        seq_base = toupper(allele_bases[a]),
        reads = depth[hit],
        is_biallelic = site_biallelic
      )
    }
    if (r > 0) {
      chunks[[i]] <- dplyr::bind_rows(rows[seq_len(r)])
      n_simple_sites <- n_simple_sites + 1L
    }
  }

  if (length(chunks) == 0L || all(vapply(chunks, is.null, logical(1)))) {
    snp_calls <- tibble::tibble(
      specimen_name = character(),
      target_name = character(),
      chrom = character(),
      pos = integer(),
      snp_name = character(),
      strand = character(),
      ref_base = character(),
      seq_base = character(),
      reads = integer(),
      is_biallelic = logical()
    )
  } else {
    snp_calls <- dplyr::bind_rows(chunks) |>
      dplyr::arrange(
        .data$specimen_name, .data$chrom, .data$pos, .data$seq_base
      )
  }

  if ("target_name" %in% colnames(snp_calls) && all(is.na(snp_calls$target_name))) {
    snp_calls <- dplyr::select(snp_calls, -"target_name")
  }

  if (!is.null(snp_calls_output)) {
    readr::write_tsv(snp_calls, snp_calls_output)
  }

  if (isTRUE(verbose)) {
    message(sprintf(
      "Wrote %d SNP calls (%d simple sites, %d samples) to %s (skipped %d non-SNP%s)",
      nrow(snp_calls), n_simple_sites, length(samples),
      if (is.null(snp_calls_output)) "<memory>" else snp_calls_output,
      n_skipped_complex,
      if (isTRUE(biallelic)) {
        sprintf(", %d multi-allelic", n_skipped_multiallelic)
      } else {
        ""
      }
    ))
  }
  snp_calls
}
