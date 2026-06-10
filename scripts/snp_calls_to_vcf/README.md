# snp\_calls\_to\_vcf.R

Contents:
* [Tool Information](#tool-information)
* [Script Usage](#script-usage)


## Tool Information

Build a VCF from the SNP pileup calls produced by `pileup_specific_snps.R`
(either the raw `snp_calls.tsv.gz` or the `collapsed_snp_calls.tsv.gz`). The VCF
carries per-sample allelic depths (`FORMAT/AD`) so downstream tools (e.g.
`moimix::getFws`) can compute within-host statistics.

Reads are summed per `(specimen, SNP, allele)` across overlapping targets (a real
sum for the raw file, a no-op for the collapsed file). `REF`/`ALT` are written on
the forward (genome) strand -- minus-strand SNPs are reverse-complemented using
the `strand` column. Multi-allelic sites are emitted with comma-separated `ALT`
and multi-value `AD` (use `--biallelic` to keep only single-ALT sites);
monomorphic sites are skipped. `GT` is derived from `AD` at `--ploidy`, and
`##contig` lengths come from the reference genome FASTA.


## Script Usage

```bash
Usage: ./snp_calls_to_vcf.R [options]


Options:
	--snp_calls=SNP_CALLS
		TSV of SNP calls from pileup_specific_snps.R (raw snp_calls.tsv.gz or collapsed_snp_calls.tsv.gz). Required columns: specimen_name, chrom, pos, snp_name, strand, ref_base, seq_base, reads

	--genome=GENOME
		Reference genome FASTA; used to write accurate ##contig=<ID=,length=> headers

	--vcf_output=VCF_OUTPUT
		Output VCF path; gzip-compressed if it ends in .gz

	--biallelic
		Keep only biallelic sites (REF + exactly one ALT); multi-allelic sites are skipped

	--ploidy=PLOIDY
		Ploidy used to render the GT field [default 2]

	--gt_min_reads=GT_MIN_READS
		Minimum reads for an allele to count as present when deriving GT [default 1]

	--overwrite
		Overwrite --vcf_output if it already exists

	--verbose
		Print a summary message when finished (otherwise run silently)

	-h, --help
		Show this help message and exit
```


```
./scripts/snp_calls_to_vcf/snp_calls_to_vcf.R --snp_calls data/example_snp_calls.tsv --genome Pf3D7.fasta --vcf_output example_snp_calls.vcf
```
