test_that("derive_gds_path replaces vcf suffixes", {
  expect_equal(derive_gds_path("x.vcf"), "x.gds")
  expect_equal(derive_gds_path("x.vcf.gz"), "x.gds")
  expect_equal(derive_gds_path("x.VCF.GZ"), "x.gds")
  expect_equal(derive_gds_path("x.bcf"), "x.bcf.gds")
  expect_equal(derive_gds_path("x.vcf", "custom.gds"), "custom.gds")
})

test_that("gds_needs_conversion respects overwrite and mtimes", {
  vcf <- tempfile(fileext = ".vcf")
  gds <- tempfile(fileext = ".gds")
  on.exit(unlink(c(vcf, gds)), add = TRUE)
  writeLines("##fileformat=VCFv4.2", vcf)
  expect_true(gds_needs_conversion(vcf, gds, overwrite = FALSE))
  expect_true(gds_needs_conversion(vcf, gds, overwrite = TRUE))
  writeLines("gds", gds)
  Sys.setFileTime(gds, file.mtime(vcf) + 10)
  expect_false(gds_needs_conversion(vcf, gds, overwrite = FALSE))
  expect_true(gds_needs_conversion(vcf, gds, overwrite = TRUE))
})

test_that("calculate_fws_from_vcf requires SeqArray and moimix", {
  err <- tryCatch(
    calculate_fws_from_vcf(vcf = tempfile()),
    error = function(e) conditionMessage(e)
  )
  expect_true(grepl("SeqArray|moimix|not found", err))
})

test_that("calculate_fws_from_vcf errors when the VCF is missing", {
  skip_if_not_installed("SeqArray")
  skip_if_not_installed("moimix")
  missing <- file.path(tempdir(), "no-such-file.vcf")
  expect_error(
    calculate_fws_from_vcf(vcf = missing),
    "not found"
  )
})

# ---- multi-allelic handling ------------------------------------------------

# Minimal VCF writer. `sites` is a list of
# list(ref = , alt = character(), ad = matrix[n_samples, 1 + length(alt)]).
write_ad_vcf <- function(path, samples, sites) {
  gt <- function(counts) {
    nz <- which(counts > 0) - 1L
    if (length(nz) == 0) {
      "./."
    } else if (length(nz) == 1) {
      paste0(nz, "/", nz)
    } else {
      paste0(nz[1], "/", nz[2])
    }
  }
  con <- file(path, "wt")
  on.exit(close(con), add = TRUE)
  writeLines(c(
    "##fileformat=VCFv4.2",
    '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">',
    '##FORMAT=<ID=AD,Number=R,Type=Integer,Description="Allelic depths">',
    '##contig=<ID=chr1,length=100000>',
    paste(c("#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO",
            "FORMAT", samples), collapse = "\t")
  ), con)
  for (i in seq_along(sites)) {
    s <- sites[[i]]
    cells <- vapply(seq_along(samples), function(j) {
      counts <- s$ad[j, ]
      paste0(gt(counts), ":", paste(counts, collapse = ","))
    }, character(1))
    writeLines(
      paste(c("chr1", i * 100L, ".", s$ref, paste(s$alt, collapse = ","),
              "100", "PASS", ".", "GT:AD", cells), collapse = "\t"),
      con
    )
  }
  invisible(path)
}

# 24 biallelic sites spanning a range of allele fractions, then 6 tri-allelic
# ones whose third allele is exactly what moimix::getFws() mishandles.
fws_test_sites <- function(n_samples) {
  sites <- list()
  for (k in 1:24) {
    f <- k / 25
    ad <- t(vapply(seq_len(n_samples), function(j) {
      alt <- as.integer(round(1000 * (((j + k) %% 5) / 4) * f))
      c(1000L - alt, alt)
    }, integer(2)))
    sites[[length(sites) + 1L]] <- list(ref = "A", alt = "C", ad = ad)
  }
  for (k in 1:6) {
    ad <- t(vapply(seq_len(n_samples), function(j) {
      switch(as.character((j + k) %% 3),
        "0" = c(0L, 0L, 1000L),   # clonal on the 2nd ALT: moimix scores Hs = 1
        "1" = c(1000L, 0L, 0L),
        c(0L, 500L, 500L)
      )
    }, integer(3)))
    sites[[length(sites) + 1L]] <- list(ref = "A", alt = c("C", "G"), ad = ad)
  }
  sites
}

test_that("calculate_fws_from_vcf drops multi-allelic sites and matches a biallelic-only callset", {
  skip_if_not_installed("SeqArray")
  skip_if_not_installed("moimix")

  samples <- paste0("S", 1:6)
  sites <- fws_test_sites(length(samples))
  is_multi <- vapply(sites, function(s) length(s$alt) > 1L, logical(1))

  dir <- tempfile("fws")
  dir.create(dir)
  on.exit(unlink(dir, recursive = TRUE), add = TRUE)
  all_vcf <- file.path(dir, "all.vcf")
  bi_vcf <- file.path(dir, "bi.vcf")
  write_ad_vcf(all_vcf, samples, sites)
  write_ad_vcf(bi_vcf, samples, sites[!is_multi])

  expect_warning(
    got <- calculate_fws_from_vcf(all_vcf, output = file.path(dir, "all.tsv")),
    sprintf("%d of %d sites are multi-allelic", sum(is_multi), length(sites))
  )
  want <- calculate_fws_from_vcf(bi_vcf, output = file.path(dir, "bi.tsv"))

  # Dropping the records inside the function must equal never having had them.
  expect_equal(got$specimen_name, want$specimen_name)
  expect_equal(got$fws, want$fws)
})

test_that("calculate_fws_from_vcf errors when no site is biallelic", {
  skip_if_not_installed("SeqArray")
  skip_if_not_installed("moimix")

  samples <- paste0("S", 1:6)
  sites <- Filter(function(s) length(s$alt) > 1L, fws_test_sites(length(samples)))
  dir <- tempfile("fws")
  dir.create(dir)
  on.exit(unlink(dir, recursive = TRUE), add = TRUE)
  vcf <- file.path(dir, "multi_only.vcf")
  write_ad_vcf(vcf, samples, sites)

  expect_error(
    calculate_fws_from_vcf(vcf, output = file.path(dir, "x.tsv")),
    "No biallelic sites"
  )
})
