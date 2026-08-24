test_that("get_aln_pos_per_real_pos skips gap characters", {
  expect_equal(get_aln_pos_per_real_pos("AC-GT", 1), 1)
  expect_equal(get_aln_pos_per_real_pos("AC-GT", 3), 4)
  expect_equal(get_aln_pos_per_real_pos("AC-GT", 4), 5)
})

test_that("add_intersected_features_to_ref_bed records contained features", {
  ref_bed <- tibble::tibble(
    `#chrom` = "chr1",
    start = 0,
    end = 10,
    target_name = "t1"
  )
  feats <- tibble::tibble(
    `#chrom` = c("chr1", "chr1", "chr2"),
    start = c(2, 9, 0),
    end = c(3, 12, 1)
  )
  out <- add_intersected_features_to_ref_bed(ref_bed, feats, "intersected")
  expect_equal(out$intersected, "1")
  covered <- add_covered_by_target_to_features(feats, ref_bed)
  expect_equal(covered$covered_by_target, c("t1", "uncovered", "uncovered"))
})

test_that("pileup_specific_snps extracts a matching SNP base", {
  skip_if_not_installed("Biostrings")
  skip_if_not_installed("pwalign")

  allele_table <- tibble::tibble(
    specimen_name = "s1",
    target_name = "t1",
    reads = 10,
    seq = "ATGAAATTT"
  )
  ref_bed <- tibble::tibble(
    `#chrom` = "chr1",
    start = 0,
    end = 9,
    target_name = "t1",
    length = 9,
    strand = "+",
    ref_seq = "ATGAAATTT"
  )
  snps <- tibble::tibble(
    `#chrom` = "chr1",
    start = 3,
    end = 4,
    name = "snpA",
    length = 1,
    strand = "+"
  )
  out_dir <- tempfile()
  on.exit(unlink(out_dir, recursive = TRUE), add = TRUE)

  result <- pileup_specific_snps(
    allele_table = allele_table,
    ref_bed = ref_bed,
    snps_of_interest = snps,
    output_directory = out_dir
  )
  expect_true(file.exists(file.path(out_dir, "snp_calls.tsv.gz")))
  expect_equal(nrow(result$snp_calls), 1)
  expect_equal(result$snp_calls$ref_base, "A")
  expect_equal(result$snp_calls$seq_base, "A")
  expect_equal(result$snp_calls$pos, 3)
})

test_that("pileup_specific_snps refuses an existing output directory", {
  skip_if_not_installed("Biostrings")
  skip_if_not_installed("pwalign")
  tmp <- tempfile()
  dir.create(tmp)
  on.exit(unlink(tmp, recursive = TRUE), add = TRUE)
  expect_error(
    pileup_specific_snps(
      allele_table = tibble::tibble(
        specimen_name = "s1",
        target_name = "t1",
        reads = 1,
        seq = "AAA"
      ),
      ref_bed = tibble::tibble(
        `#chrom` = "chr1",
        start = 0,
        end = 3,
        target_name = "t1",
        length = 3,
        strand = "+",
        ref_seq = "AAA"
      ),
      snps_of_interest = tibble::tibble(
        `#chrom` = "chr1",
        start = 0,
        end = 1,
        name = "s",
        length = 1,
        strand = "+"
      ),
      output_directory = tmp
    ),
    "already exist"
  )
})
