Rscript scripts/slaf_from_stave_mlaf/slaf_from_stave_mlaf.R --mlaf_input data/example_mlaf.tsv --output single_locus_allele_freq.tsv
Rscript scripts/filter_biallelic_calls/filter_biallelic_calls.R --amino_acid_calls data/example2_amino_acid_calls.tsv --out biallelic.tsv
Rscript scripts/allele_per_locus_summary/allele_per_locus_summary.R --allele_table data/example2_allele_table.tsv
Rscript scripts/add_ref_seq_to_ref_bed_table/add_ref_seqs_with_full_genome_ref_fasta.R --ref_bed data/example_PMO_insert_locs_of_panel.bed --genome_fasta ~/Downloads/PkPfPmPoPv.fasta --out example.bed 
Rscript scripts/add_ref_seq_to_ref_bed_table/add_ref_seqs_with_targeted_ref_fasta.R --ref_bed data/example_PMO_insert_locs_of_panel.bed --target_fasta data/example_PMO_insert_locs_of_panel_refseqs.fasta --out test_out.bed 
Rscript scripts/translate_loci_of_interest/translate_loci_of_interest.R --allele_table data/example2_allele_table.tsv --ref_bed data/example_PMO_insert_locs_of_panel.bed --loci_of_interest data/example_principal_resistance_marker_info_table.bed --output_directory tmp_dir --overwrite_dir
Rscript scripts/THEREALMcCOIL_wrapper/THEREALMcCOIL_wrapper.R --model categorical --snp_calls_input data/example_snp_calls.tsv --slaf_output test_data_out_slaf.tsv --coi_output test_data_output_coi.tsv 
Rscript scripts/calc_slaf_based_on_mhap_freqs/slaf_from_mhaps_freqs.R --mhaps_slaf_fnp data/example_mhaps_slaf.tsv --loci_of_interest_per_microhaps_fnp data/example_loci_of_interest_for_target_for_microhap.tsv --slaf_output slaf.tsv 
Rscript scripts/snpslice_wrapper/snpslice_wrapper.R --allele_table data/example_amino_acid_calls.tsv --loci_groups_input data/example_loci_groups.tsv --mlaf_output mlaf.tsv --coi_output coi.tsv
Rscript scripts/pileup_specific_snps/pileup_specific_snps.R --allele_table data/example2_allele_table.tsv --ref_bed data/example_PMO_insert_locs_of_panel.bed --snps_of_interest data/MAD4HATTER_coveredSnps.bed --output_directory tmp_output --overwrite_dir
Rscript scripts/per_locus_popgen_summary_wrapper/per_locus_tajima_d_summary_wrapper.R --allele_table data/example2_allele_table.tsv
Rscript scripts/multilocus_prevfreq_naive/multilocus_prevfreq_naive.R --aa_table data/example_amino_acid_calls.tsv --loci_groups_input data/example_loci_groups.tsv --output_path mlafp.tsv 
Rscript scripts/moire_wrapper/moire_wrapper.R --allele_table data/example2_allele_table.tsv --mcmc_results_output moire_output.tsv
Rscript scripts/MultiLociBiallelicModel_wrapper/MultiLociBiallelicModel_wrapper.R --aa_calls data/example_amino_acid_calls.tsv --loci_group_table data/example_loci_groups.tsv --mlaf_output mlaf.tsv
Rscript scripts/malariaem_wrapper/malariaem_wrapper.R  --allele_table data/example2_allele_table.tsv --subset_targets TRUE --target_groups data/example_target_groups.tsv --freq_output freq_out.tsv --phase_out phase_out.tsv
Rscript scripts/IDM_wrapper/IDM_wrapper.R --allele_table_input data/example_allele_table.tsv --model IDM --slaf_output out.tsv 
Rscript scripts/FreqEstimationModel_wrapper/FreqEstimationModel_wrapper.R --aa_calls data/example_amino_acid_calls.tsv --coi data/example_coi_table.tsv --groups data/example_loci_groups.tsv --mlaf_output output.tsv
Rscript scripts/filter_to_highest_diversity_independent_snp_call/filter_to_highest_diversity_independent_snp_call.R --snp_table_in data/example_collapsed_snp_calls.tsv --snp_table_out tmp.tsv.gz 
Rscript scripts/estimate_coi_naive/estimate_coi_naive.R --input_path data/example_allele_table.tsv --output_path coi_table.tsv --method integer_method --integer_threshold 5
Rscript scripts/estimate_allele_prevalence_naive/estimate_allele_prevalence_naive.R --aa_calls data/example_amino_acid_calls.tsv --output prevalence.tsv
Rscript scripts/estimate_allele_frequency_naive/estimate_allele_frequency_naive.R --aa_calls data/example_amino_acid_calls.tsv --method read_count_prop --output allele_freqs.tsv
Rscript scripts/dcifer_slaf_wrapper/dcifer_slaf_wrapper.R --allele_table data/example2_amino_acid_calls.tsv --target_name_col aa_locus --target_value_col aa --slaf_output slaf.tsv
Rscript scripts/dcifer_ibd_wrapper/dcifer_ibd_wrapper.R --allele_table data/example2_allele_table.tsv --allele_freq_table data/example_slaf_mhap.tsv --btwn_host_rel_output btwn_host_rel.tsv
Rscript scripts/count_samples_by_coi/count_samples_by_coi.R --coi_calls data/example_coi_table.tsv --output coi_counts.tsv
Rscript scripts/coiaf_wrapper/coiaf_wrapper.R --snp_data data/example_collapsed_snp_calls.tsv --output coif_output
Rscript scripts/allele_per_locus_summary/allele_per_locus_summary.R --allele_table data/example_allele_table.tsv
