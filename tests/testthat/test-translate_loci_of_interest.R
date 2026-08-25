test_that("translate_loci_of_interest translates a plus-strand codon", {
  skip_if_not_installed("Biostrings")
  skip_if_not_installed("pwalign")

  allele_table <- tibble::tibble(
    specimen_name = "s1",
    target_name = "t1",
    reads = 5,
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
  loci <- tibble::tibble(
    `#chrom` = "chr1",
    start = 0,
    end = 3,
    name = "aa1",
    length = 3,
    strand = "+",
    gene = "g",
    gene_id = "G1",
    aa_position = 1
  )
  out_dir <- tempfile()
  on.exit(unlink(out_dir, recursive = TRUE), add = TRUE)

  result <- translate_loci_of_interest(
    allele_table = allele_table,
    ref_bed = ref_bed,
    loci_of_interest = loci,
    output_dir = out_dir
  )
  expect_true(
    file.exists(file.path(out_dir, "amino_acid_calls.tsv.gz"))
  )
  expect_equal(result$amino_acid_calls$aa, "M")
  expect_equal(result$amino_acid_calls$ref_aa, "M")
  expect_equal(result$amino_acid_calls$codon, "ATG")
  expect_equal(result$amino_acid_calls$aa_locus, "G1:1")
})

test_that("translate_loci_of_interest rejects loci that are not length 3", {
  skip_if_not_installed("Biostrings")
  skip_if_not_installed("pwalign")
  out_dir <- tempfile()
  on.exit(unlink(out_dir, recursive = TRUE), add = TRUE)
  expect_error(
    translate_loci_of_interest(
      allele_table = tibble::tibble(
        specimen_name = "s1",
        target_name = "t1",
        reads = 1,
        seq = "ATG"
      ),
      ref_bed = tibble::tibble(
        `#chrom` = "chr1",
        start = 0,
        end = 3,
        target_name = "t1",
        length = 3,
        strand = "+",
        ref_seq = "ATG"
      ),
      loci_of_interest = tibble::tibble(
        `#chrom` = "chr1",
        start = 0,
        end = 1,
        name = "bad",
        length = 1,
        strand = "+",
        gene = "g",
        gene_id = "G1",
        aa_position = 1
      ),
      output_dir = out_dir
    ),
    "length 3"
  )
})
