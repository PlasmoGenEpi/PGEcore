# Package index

## Getting started helpers

Lightweight summaries and COI utilities with no specialised Suggests.

- [`count_samples_by_coi()`](https://plasmogenepi.github.io/PGEcore/reference/count_samples_by_coi.md)
  : Count specimens by complexity of infection (COI)
- [`estimate_coi_naive()`](https://plasmogenepi.github.io/PGEcore/reference/estimate_coi_naive.md)
  : Estimate COI using naive allele-count methods
- [`estimate_allele_frequency_naive()`](https://plasmogenepi.github.io/PGEcore/reference/estimate_allele_frequency_naive.md)
  : Estimate allele frequency naively from AA or microhaplotype calls
- [`estimate_allele_prevalence_naive()`](https://plasmogenepi.github.io/PGEcore/reference/estimate_allele_prevalence_naive.md)
  : Estimate allele prevalence naively from AA or microhaplotype calls
- [`allele_per_locus_summary()`](https://plasmogenepi.github.io/PGEcore/reference/allele_per_locus_summary.md)
  : Summarize alleles per locus from an allele table

## Format converters and filters

- [`filter_biallelic_calls()`](https://plasmogenepi.github.io/PGEcore/reference/filter_biallelic_calls.md)
  : Filter amino acid calls to biallelic loci
- [`filter_to_highest_diversity_independent_snp_call()`](https://plasmogenepi.github.io/PGEcore/reference/filter_to_highest_diversity_independent_snp_call.md)
  : Filter SNPs to highest-diversity loci spaced by a minimum distance
- [`slaf_from_mhaps_freqs()`](https://plasmogenepi.github.io/PGEcore/reference/slaf_from_mhaps_freqs.md)
  : Calculate single-locus allele frequencies from microhaplotype
  frequencies
- [`slaf_from_stave_mlaf()`](https://plasmogenepi.github.io/PGEcore/reference/slaf_from_stave_mlaf.md)
  : Convert STAVE multi-locus allele frequencies to single-locus
  frequencies
- [`convert_single_locus_table_to_stave()`](https://plasmogenepi.github.io/PGEcore/reference/convert_single_locus_table_to_stave.md)
  : Convert a single-locus table to STAVE-style variant identifiers
- [`multilocus_prevfreq_naive()`](https://plasmogenepi.github.io/PGEcore/reference/multilocus_prevfreq_naive.md)
  : Estimate multilocus prevalence and frequency with naive phasing
- [`multilocus_prevfreq_naive_variantstring()`](https://plasmogenepi.github.io/PGEcore/reference/multilocus_prevfreq_naive_variantstring.md)
  : Estimate multilocus prevalence and frequency with variantstring
- [`snp_calls_to_vcf()`](https://plasmogenepi.github.io/PGEcore/reference/snp_calls_to_vcf.md)
  : Build a VCF from pileup SNP calls
- [`vcf_to_snp_calls()`](https://plasmogenepi.github.io/PGEcore/reference/vcf_to_snp_calls.md)
  : Convert a VCF with FORMAT/AD into pileup-style SNP calls
- [`add_ref_seqs_with_targeted_ref_fasta()`](https://plasmogenepi.github.io/PGEcore/reference/add_ref_seqs_with_targeted_ref_fasta.md)
  : Add reference sequences from a targeted FASTA onto a panel BED table
- [`add_ref_seqs_with_full_genome_ref_fasta()`](https://plasmogenepi.github.io/PGEcore/reference/add_ref_seqs_with_full_genome_ref_fasta.md)
  : Add reference sequences extracted from a genome FASTA onto a panel
  BED table

## Sequence and panel tools

Require Biostrings and related Suggests (and often MSA binaries on
PATH).

- [`pileup_specific_snps()`](https://plasmogenepi.github.io/PGEcore/reference/pileup_specific_snps.md)
  : Pile up specific SNPs covered by microhaplotype sequences
- [`translate_loci_of_interest()`](https://plasmogenepi.github.io/PGEcore/reference/translate_loci_of_interest.md)
  : Translate loci of interest from microhaplotype sequences
- [`per_locus_popgen_summary()`](https://plasmogenepi.github.io/PGEcore/reference/per_locus_popgen_summary.md)
  : Per-locus nucleotide diversity, segregating sites, and Tajima's D

## Specialist wrappers

Optional packages must be installed separately (see README Suggests
table). Prefer these file/CLI entry points.

- [`coiaf_wrapper()`](https://plasmogenepi.github.io/PGEcore/reference/coiaf_wrapper.md)
  : Estimate COI with coiaf from SNP-call and output paths
- [`dcifer_slaf_wrapper()`](https://plasmogenepi.github.io/PGEcore/reference/dcifer_slaf_wrapper.md)
  : Estimate single-locus allele frequencies with Dcifer
- [`dcifer_ibd_wrapper()`](https://plasmogenepi.github.io/PGEcore/reference/dcifer_ibd_wrapper.md)
  : Estimate IBD-based relatedness with Dcifer
- [`moire_wrapper()`](https://plasmogenepi.github.io/PGEcore/reference/moire_wrapper.md)
  : Run MOIRe from allele-table and output paths
- [`malariaem_wrapper()`](https://plasmogenepi.github.io/PGEcore/reference/malariaem_wrapper.md)
  : Run malaria.em from allele-table and output paths
- [`snpslice_wrapper()`](https://plasmogenepi.github.io/PGEcore/reference/snpslice_wrapper.md)
  : Estimate multilocus allele frequency and COI with SNP-Slice
- [`FreqEstimationModel_wrapper()`](https://plasmogenepi.github.io/PGEcore/reference/FreqEstimationModel_wrapper.md)
  : Estimate multilocus allele frequencies with FreqEstimationModel
- [`calculate_fws_from_vcf()`](https://plasmogenepi.github.io/PGEcore/reference/calculate_fws_from_vcf.md)
  : Calculate within-host Fws from a VCF via moimix
- [`IDM_wrapper()`](https://plasmogenepi.github.io/PGEcore/reference/IDM_wrapper.md)
  : Estimate single-locus allele frequencies with the Incomplete Data
  Model
- [`MultiLociBiallelicModel_wrapper()`](https://plasmogenepi.github.io/PGEcore/reference/MultiLociBiallelicModel_wrapper.md)
  : Estimate multilocus haplotype frequencies with
  MultiLociBiallelicModel
- [`THEREALMcCOIL_wrapper()`](https://plasmogenepi.github.io/PGEcore/reference/THEREALMcCOIL_wrapper.md)
  : Estimate COI and allele frequencies with THEREALMcCOIL

## In-memory helpers

Optional `run_*` APIs for a few wrappers when data are already in R.

- [`run_coiaf()`](https://plasmogenepi.github.io/PGEcore/reference/run_coiaf.md)
  : Estimate complexity of infection (COI) using coiaf
- [`run_moire()`](https://plasmogenepi.github.io/PGEcore/reference/run_moire.md)
  : Run MOIRe MCMC analysis
- [`run_malariaem()`](https://plasmogenepi.github.io/PGEcore/reference/run_malariaem.md)
  : Run malaria.em and write frequency and phase summaries
