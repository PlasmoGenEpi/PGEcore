test_that("add_ref_seqs_with_full_genome_ref_fasta checks Biostrings and duplicates", {
  err <- tryCatch(
    add_ref_seqs_with_full_genome_ref_fasta(
      ref_bed = data.frame(x = 1),
      genome_fasta = tempfile()
    ),
    error = function(e) conditionMessage(e)
  )
  expect_true(grepl("Biostrings|Missing required columns", err))

  bed <- tibble::tibble(
    `#chrom` = c("chr1", "chr1"),
    start = c(0, 1),
    end = c(2, 3),
    target_name = c("t1", "t1"),
    length = c(2, 2),
    strand = c("+", "+")
  )
  skip_if_not_installed("Biostrings")
  expect_error(
    add_ref_seqs_with_full_genome_ref_fasta(bed, tempfile()),
    "multiple times: t1"
  )
})

test_that("add_ref_seqs_with_full_genome_ref_fasta extracts plus and minus strand intervals", {
  skip_if_not_installed("Biostrings")

  fasta <- tempfile(fileext = ".fa")
  on.exit(unlink(fasta), add = TRUE)
  writeLines(c(">chr1 comment", "ACGTACGT"), fasta)

  bed <- tibble::tibble(
    `#chrom` = c("chr1", "chr1"),
    start = c(1, 1),
    end = c(5, 5),
    target_name = c("plus", "minus"),
    length = c(4, 4),
    strand = c("+", "-")
  )

  out <- add_ref_seqs_with_full_genome_ref_fasta(bed, fasta)
  expect_equal(out$ref_seq[out$target_name == "plus"], "CGTA")
  expect_equal(out$ref_seq[out$target_name == "minus"], "TACG")

  missing_chrom <- bed[1, ]
  missing_chrom[["#chrom"]] <- "chrX"
  expect_error(
    add_ref_seqs_with_full_genome_ref_fasta(missing_chrom, fasta),
    "not in"
  )
})
