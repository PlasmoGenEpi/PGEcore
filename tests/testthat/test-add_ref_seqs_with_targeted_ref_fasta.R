test_that("validate_ref_bed_table requires columns and types", {
  expect_error(
    validate_ref_bed_table(data.frame(target_name = "t"), "ref_bed"),
    "Missing required columns"
  )

  ok <- tibble::tibble(
    `#chrom` = "chr1",
    start = 0,
    end = 4,
    target_name = "t1",
    length = 4,
    strand = "+"
  )
  expect_true(validate_ref_bed_table(ok))
})

test_that("stop_on_duplicate_names reports repeated target_name values", {
  expect_error(
    stop_on_duplicate_names(c("a", "a", "b"), "demo.bed"),
    "multiple times: a"
  )
  expect_true(stop_on_duplicate_names(c("a", "b"), "demo.bed"))
})

test_that("add_ref_seqs_with_targeted_ref_fasta checks Biostrings and overwrite", {
  err <- tryCatch(
    add_ref_seqs_with_targeted_ref_fasta(
      ref_bed = data.frame(x = 1),
      target_fasta = tempfile()
    ),
    error = function(e) conditionMessage(e)
  )
  expect_true(grepl("Biostrings|Missing required columns", err))

  tmp <- tempfile()
  on.exit(unlink(tmp), add = TRUE)
  writeLines("x", tmp)
  expect_error(
    add_ref_seqs_with_targeted_ref_fasta(
      ref_bed = data.frame(x = 1),
      target_fasta = tempfile(),
      output = tmp,
      overwrite = FALSE
    ),
    "already exists"
  )
})

test_that("add_ref_seqs_with_targeted_ref_fasta joins FASTA records by target_name", {
  skip_if_not_installed("Biostrings")

  fasta <- tempfile(fileext = ".fa")
  on.exit(unlink(fasta), add = TRUE)
  writeLines(c(">t1", "ACGT", ">t2", "GG"), fasta)

  bed <- tibble::tibble(
    `#chrom` = c("chr1", "chr1"),
    start = c(0, 10),
    end = c(4, 12),
    target_name = c("t1", "t2"),
    length = c(4, 2),
    strand = c("+", "+"),
    ref_seq = c("old", "old")
  )

  out <- add_ref_seqs_with_targeted_ref_fasta(bed, fasta)
  expect_equal(out$ref_seq, c("ACGT", "GG"))

  missing <- bed[1, ]
  missing$target_name <- "missing_locus"
  expect_error(
    add_ref_seqs_with_targeted_ref_fasta(missing, fasta),
    "missing from the fasta file"
  )
})
