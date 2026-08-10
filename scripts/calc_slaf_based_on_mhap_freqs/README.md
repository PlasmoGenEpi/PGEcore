# Calculate Allele Frequency Associated with Microhaplotype based on Microhaplotype Frequencies

TODO: come back to this README

Contents:

-   [Tool Information](#tool-information)
-   [Script Usage](#script-usage)

## Tool Information {#tool-information}

Taking amino acid translated loci associated with microhaplotype sequences and calculate their allele frequencies based on the allele frequencies of the microhaplotypes.

## Script Usage {#script-usage}

Script requires two inputs (the translated loci associated with microhaplotypes and the allele frequencies for the microhaplotypes) and can write out the frequencies collapsed across any overlapping targets (and optionally can also export per target_name as well)

```         
Usage: ./slaf_from_mhaps_freqs.R [options]

Options:
    --mhaps_slaf_fnp=MHAPS_SLAF_FNP
        TSV containing the columns: target_name, seq, freq, sample_total. The target_name and seq columns should match up with the columns in loci_of_interest_per_microhaps_fnp

    --loci_of_interest_per_microhaps_fnp=LOCI_OF_INTEREST_PER_MICROHAPS_FNP
        TSV containing the columns: target_name, seq, gene_id, aa_position, aa. The target_name and seq columns should match up with the columns in mhaps_slaf_fnp

    --slaf_output=SLAF_OUTPUT
        the output for the single locus allele frequency, will collapse frequencies across overlapping targets

    --per_target_slaf_output=PER_TARGET_SLAF_OUTPUT
        optional output for the single locus allele frequency calculated per 
        target

    -h, --help
        Show this help message and exit
```

### Examples

The two inputs can be calculated using other scripts in PGEcore

```         
scripts/dcifer_slaf_wrapper/dcifer_slaf_wrapper.R --allele_table \
    data/example2_allele_table.tsv \
    --slaf_output mhaps_slaf.tsv

./scripts/translate_loci_of_interest/translate_loci_of_interest.R \
    --output_directory translate_output --allele_table data/example2_allele_table.tsv  \
    --ref_bed data/example_PMO_insert_locs_of_panel.bed  \
    --loci_of_interest data/example_principal_resistance_marker_info_table.bed \
    --overwrite_dir

./scripts/calc_slaf_based_on_mhap_freqs/slaf_from_mhaps_freqs.R --mhaps_slaf_fnp mhaps_slaf.tsv \ 
    --loci_of_interest_per_microhaps_fnp translate_output/loci_of_interest_for_target_for_microhap.tsv.gz \
    --slaf_output slaf.tsv
```

Allele frequencies can be calculated in any way to be given to this script

```         
./scripts/calc_slaf_based_on_mhap_freqs/slaf_from_mhaps_freqs.R --mhaps_slaf_fnp data/example_mhaps_slaf.tsv \ 
    --loci_of_interest_per_microhaps_fnp data/example_loci_of_interest_for_target_for_microhap.tsv \
    --slaf_output slaf.tsv
```
