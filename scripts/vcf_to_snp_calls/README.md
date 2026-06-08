# vcf\_to\_snp\_calls.R

Contents:
* [Tool Information](#tool-information)
* [Script Usage](#script-usage)


## Tool Information

The reverse of `snp_calls_to_vcf.R`: read a VCF and emit per-sample SNP calls in
a format as close as possible to `pileup_specific_snps.R`'s `snp_calls.tsv.gz`.
Read counts are sourced from `FORMAT/AD` (allelic depths), so the input VCF must
carry `AD` -- e.g. one built by `snp_calls_to_vcf.R`.

One row is written per `(specimen, SNP, observed allele)` for every allele whose
`AD` is at least `--min_reads` (default 1). Only **simple** SNP alleles are kept
(single-base `REF` and `ALT` in `{A,C,G,T}`); MNPs, indels, and symbolic alleles
(e.g. `*`, `<NON_REF>`) are skipped. With `--biallelic`, only sites with exactly
one simple `ALT` (REF + 1 ALT) are kept.

The emitted columns mirror the pileup SNP-calls shape: `specimen_name`,
`chrom`, `pos` (0-based), `snp_name`, `strand`, `ref_base`, `seq_base`, `reads`,
`is_biallelic`. Notes on the round-trip:

* `pos` is written 0-based (VCF `POS` minus 1), matching the pileup output.
* `snp_name` uses the VCF `ID`; when `ID` is `.` it falls back to
  `chrom-start-end` (0-based start), matching the pileup naming.
* `strand` is always `+` -- `REF`/`ALT` are on the forward (genome) strand.
* `target_name` is emitted **only** when the variant's `INFO` carries a
  `TARGET=` key; otherwise the column is omitted entirely (a plain VCF has no
  microhaplotype/target context).
* Columns with no source in a VCF (`seq`, `he`, ...) are not emitted.


## Script Usage

```bash
Usage: ./vcf_to_snp_calls.R [options]


Options:
	--vcf=VCF
		Input VCF (.vcf or .vcf.gz). Must carry FORMAT/AD (allelic depths).

	--snp_calls_output=SNP_CALLS_OUTPUT
		Output SNP-calls TSV path; gzip-compressed if it ends in .gz

	--biallelic
		Keep only biallelic sites (REF + exactly one simple ALT); multi-allelic sites are skipped

	--min_reads=MIN_READS
		Minimum AD reads for an allele to be emitted as an observed call [default 1]

	--overwrite
		Overwrite --snp_calls_output if it already exists

	--verbose
		Print a summary message when finished (otherwise run silently)

	-h, --help
		Show this help message and exit
```


```
./scripts/vcf_to_snp_calls/vcf_to_snp_calls.R --vcf data/example_snp_calls.vcf --snp_calls_output example_snp_calls_from_vcf.tsv
```
