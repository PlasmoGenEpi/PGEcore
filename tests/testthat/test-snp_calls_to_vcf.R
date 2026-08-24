test_that("revcomp reverse-complements DNA strings", {
  expect_equal(revcomp("A"), "T")
  expect_equal(revcomp("AC"), "GT")
  expect_equal(revcomp(c("A", "Nn")), c("T", "nN"))
})

test_that("derive_gt pads, truncates, and marks missing alleles", {
  expect_equal(derive_gt(c(10, 5), ploidy = 2L, min_reads = 1L), "0/1")
  expect_equal(derive_gt(c(10, 0), ploidy = 2L, min_reads = 1L), "0/0")
  expect_equal(derive_gt(c(0, 0), ploidy = 2L, min_reads = 1L), "./.")
  expect_equal(derive_gt(c(1, 10, 5), ploidy = 2L, min_reads = 1L), "1/2")
  expect_equal(derive_gt(c(3, 3), ploidy = 2L, min_reads = 10L), "./.")
})

test_that("snp_calls_to_vcf validates columns and Biostrings", {
  err <- tryCatch(
    snp_calls_to_vcf(
      snp_calls = data.frame(x = 1),
      genome = tempfile(),
      vcf_output = tempfile(fileext = ".vcf")
    ),
    error = function(e) conditionMessage(e)
  )
  expect_true(grepl("Biostrings|Missing required columns", err))
})

test_that("snp_calls_to_vcf writes a biallelic VCF when Biostrings is present", {
  skip_if_not_installed("Biostrings")

  fasta <- tempfile(fileext = ".fa")
  vcf <- tempfile(fileext = ".vcf")
  on.exit(unlink(c(fasta, vcf)), add = TRUE)
  writeLines(c(">chr1 extra description", "ACGTACGT"), fasta)

  calls <- tibble::tibble(
    specimen_name = c("s1", "s1", "s1"),
    chrom = "chr1",
    pos = 0L,
    snp_name = "chr1-0-1",
    strand = "+",
    ref_base = "A",
    seq_base = c("A", "G", "A"),
    reads = c(8L, 2L, 1L)
  )

  expect_true(
    snp_calls_to_vcf(
      snp_calls = calls,
      genome = fasta,
      vcf_output = vcf,
      overwrite = TRUE
    )
  )
  lines <- readLines(vcf)
  expect_true(any(grepl("^##fileformat=VCFv4.2", lines)))
  expect_true(any(grepl("##contig=<ID=chr1,length=8>", lines)))
  body <- lines[!startsWith(lines, "##") & !startsWith(lines, "#CHROM")]
  expect_equal(length(body), 1L)
  fields <- strsplit(body, "\t")[[1]]
  expect_equal(fields[1:5], c("chr1", "1", "chr1-0-1", "A", "G"))
  expect_match(fields[10], "^0/1:9,2:11$")
})

test_that("snp_calls_to_vcf reverse-complements minus-strand alleles", {
  skip_if_not_installed("Biostrings")

  fasta <- tempfile(fileext = ".fa")
  vcf <- tempfile(fileext = ".vcf")
  on.exit(unlink(c(fasta, vcf)), add = TRUE)
  writeLines(c(">chr1", "ACGTACGT"), fasta)

  calls <- tibble::tibble(
    specimen_name = c("s1", "s1"),
    chrom = "chr1",
    pos = 1L,
    snp_name = "snp_minus",
    strand = "-",
    ref_base = "T",
    seq_base = c("T", "C"),
    reads = c(4L, 3L)
  )

  snp_calls_to_vcf(calls, fasta, vcf, overwrite = TRUE)
  body <- readLines(vcf)
  rec <- body[!startsWith(body, "#")]
  fields <- strsplit(rec, "\t")[[1]]
  expect_equal(fields[4], "A")
  expect_equal(fields[5], "G")
})

test_that("snp_calls_to_vcf refuses to overwrite", {
  tmp <- tempfile(fileext = ".vcf")
  on.exit(unlink(tmp), add = TRUE)
  writeLines("x", tmp)
  expect_error(
    snp_calls_to_vcf(
      snp_calls = data.frame(specimen_name = "s"),
      genome = tempfile(),
      vcf_output = tmp,
      overwrite = FALSE
    ),
    "already exists"
  )
})
