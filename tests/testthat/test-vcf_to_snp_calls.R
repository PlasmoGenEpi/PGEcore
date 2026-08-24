test_that("is_simple_base accepts single DNA letters", {
  expect_true(is_simple_base("A"))
  expect_true(is_simple_base("g"))
  expect_false(is_simple_base("AT"))
  expect_false(is_simple_base("*"))
  expect_false(is_simple_base("<NON_REF>"))
})

test_that("info_target parses TARGET= from INFO", {
  expect_equal(info_target("NS=1;TARGET=amplicon1;DP=10"), "amplicon1")
  expect_equal(info_target("TARGET=only"), "only")
  expect_true(is.na(info_target("NS=1;DP=10")))
})

test_that("vcf_to_snp_calls emits 0-based pileup-style rows from FORMAT/AD", {
  vcf <- tempfile(fileext = ".vcf")
  on.exit(unlink(vcf), add = TRUE)
  writeLines(
    c(
      "##fileformat=VCFv4.2",
      "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\ts1\ts2",
      "chr1\t1\t.\tA\tG\t.\t.\tNS=2\tGT:AD:DP\t0/1:8,2:10\t0/0:5,0:5",
      "chr1\t2\tindel1\tAT\tA\t.\t.\t.\tGT:AD:DP\t0/1:1,1:2\t.:.:."
    ),
    vcf
  )

  out <- vcf_to_snp_calls(vcf)
  expect_false("target_name" %in% names(out))
  expect_equal(unique(out$strand), "+")
  expect_equal(sort(unique(out$pos)), 0L)
  expect_equal(unique(out$snp_name), "chr1-0-1")
  expect_equal(nrow(out), 3L)
  s1_alt <- out[out$specimen_name == "s1" & out$seq_base == "G", ]
  expect_equal(s1_alt$reads, 2L)
})

test_that("vcf_to_snp_calls keeps TARGET= as target_name and skips multi-allelic when requested", {
  vcf <- tempfile(fileext = ".vcf")
  on.exit(unlink(vcf), add = TRUE)
  writeLines(
    c(
      "##fileformat=VCFv4.2",
      "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\ts1",
      "chr1\t10\tsnpA\tC\tT,G\t.\t.\tTARGET=t1\tGT:AD:DP\t0/1/2:1,2,3:6",
      "chr1\t11\tsnpB\tA\tT\t.\t.\tTARGET=t2\tGT:AD:DP\t0/1:4,1:5"
    ),
    vcf
  )

  all_sites <- vcf_to_snp_calls(vcf, min_reads = 1L)
  expect_true("target_name" %in% names(all_sites))
  expect_equal(sort(unique(all_sites$snp_name)), c("snpA", "snpB"))

  biallelic <- vcf_to_snp_calls(vcf, biallelic = TRUE)
  expect_equal(unique(biallelic$snp_name), "snpB")
  expect_equal(unique(biallelic$target_name), "t2")
})

test_that("vcf_to_snp_calls requires #CHROM, AD, and overwrite protection", {
  bad <- tempfile(fileext = ".vcf")
  on.exit(unlink(bad), add = TRUE)
  writeLines(c("##fileformat=VCFv4.2", "chr1\t1\t.\tA\tG\t.\t.\t.\tGT\t0/1"), bad)
  expect_error(vcf_to_snp_calls(bad), "no #CHROM header")

  vcf <- tempfile(fileext = ".vcf")
  on.exit(unlink(vcf), add = TRUE)
  writeLines(
    c(
      "##fileformat=VCFv4.2",
      "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\ts1",
      "chr1\t1\t.\tA\tG\t.\t.\t.\tGT:DP\t0/1:10"
    ),
    vcf
  )
  expect_error(vcf_to_snp_calls(vcf), "FORMAT has no AD field")

  existing <- tempfile(fileext = ".tsv")
  on.exit(unlink(existing), add = TRUE)
  writeLines("x", existing)
  expect_error(
    vcf_to_snp_calls(vcf, snp_calls_output = existing, overwrite = FALSE),
    "already exists"
  )
})
