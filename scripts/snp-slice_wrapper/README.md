# SNP-slice

Contents: 
* [Tool Information](#tool-information)
* [Script Usage](#script-usage)

## Tool Information

The publication associated with the SNP-slice algorithm and analysis scripts can be found here [Ju et al 2024 Bioinformatics](https://academic.oup.com/bioinformatics/article/40/6/btae344/7695237?login=false) and the original code can be found [here](https://github.com/nianqiaoju/snp-slice/tree/main).

More information about this tool can be found [here](https://mrc-ide.github.io/PGEforge/tutorials/SNP-slice/background.html).

### Purpose

This package is a wrapper for SNP-Slice, a programme to phase mutlilocus genotypes. In this wrapper we use the phased assignments to population-level multi-locus allele frequencies. For example, dhfr-dhps quintuple mutation. 

## Script Usage 

### Dependencies
This package currently depends on R and the following R packages:
 - dplyr
 - tidyr
 - readr
 - optparse
 - stringr
 - tibble

SNP-slice is not an installable R package, therefore the GitHub repo must be cloned to run this code.
```
git clone https://github.com/nianqiaoju/snp-slice.git
```
### Usage 
Basic usage 
```
Rscript scripts/snp-slice_wrapper/snp-slice_wrapper.R --aa_calls \
    data/example_amino_acid_calls.tsv --snp_slice_dir path/to/snp-slice/ \
    --snpslicemain_path scripts/snp-slice_wrapper/adapted_snpslicemain.R \
    --loci_group_table data/example_loci_groups.tsv  --model 3 --gap 50 \ 
    --output output.tsv
```

#### Required Arguments

- `--aa_calls`: Path to amino acid calls (TSV format)
- `--snp_slice_dir`: Path to local copy of SNP-slice GitHub repository 
- `--snpslicemain_path`: Path to modified version of `snpslicemain.R`
- `--loci_group_table`: Path to locus groups file (TSV format)
- `--model`: Model to be used with SNP-slice. See SNP-slice documentation (int)
- `--gap`: How many iterations of no improvement before stopping
- `--output`: Path to output file 

## Input File Formats

### Amino Acid Calls (Required)

A tab-separated values (TSV) file with the following columns:

| Column | Type | Description |
|--------|------|-------------|
| `specimen_id` | String | Unique identifier for each specimen |
| `gene_id` | String | Gene identifier |
| `aa_position` | String | Amino acid position |
| `ref_aa` | String | Reference amino acid |
| `aa` | String | Amino acid call |
| `read_count` | Integer | Number of reads supporting this allele |

**Example Amino Acid data:**
```tsv
specimen_id	read_count	gene_id	aa_position	ref_aa	aa
specimen1	5	PF3D7_0417200	51	N	N
specimen1	5	PF3D7_0417200	59	C	R
specimen1	5	PF3D7_0417200	108	S	N
specimen1	5	PF3D7_0810800	437	A	A
...
```

### Locus Group Table

A tab-separated values (TSV) file to define the groups of loci to be run through SNP-slice individually. Include the following columns:

| Column | Type | Description |
|--------|------|-------------|
| `group_id` | String | Unique identifier for each group to run through |
| `gene_id` | String | Gene identifier |
| `aa_position` | Numeric | Amino acid position |

These must match what is in the amino acid file. 

**Example group data:**
```tsv
group_id	gene_id	aa_position
pfdhfr_pfdhps	PF3D7_0417200	51
pfdhfr_pfdhps	PF3D7_0417200	59
pfdhfr_pfdhps	PF3D7_0417200	108
pfdhfr_pfdhps	PF3D7_0810800	437
pfdhfr_pfdhps	PF3D7_0810800	540
pfdhfr	PF3D7_0417200	51
pfdhfr	PF3D7_0417200	59
pfdhfr	PF3D7_0417200	108
pfdhps	PF3D7_0810800	437
pfdhps	PF3D7_0810800	540
```

## Output Format

The script outputs a tab-separated values (TSV) file with the following columns:

| Column | Type | Description |
|--------|------|-------------|
| `stave_string` | String | Unique identifier for the genotype |
| `frequency` | Numeric | Estimate of frequency |
| `group_id` | Numeric | The group identifier |

**Example output:**
```tsv
stave_string	frequency	group_id
51_59_108:NRS;437_540:GE	0.333333333333333	pfdhfr_pfdhps
51_59_108:IRS;437_540:GK	0.333333333333333	pfdhfr_pfdhps
51_59_108:NRN;437_540:AE	0.333333333333333	pfdhfr_pfdhps
51_59_108:NRS	0.2	pfdhfr
51_59_108:IRS	0.4	pfdhfr
51_59_108:NRN	0.4	pfdhfr
437_540:AE	0.2	pfdhps
437_540:GK	0.4	pfdhps
437_540:GE	0.4	pfdhps
```

## Functions
 - create_SNP_slice_input_df: converts amino acid tables into SNP-Slice input
 - create_ref_alt_df: converts standardized input tables into one table per
 amino acid position, with a single reference allele and a single alternate
 allele associated with each position.
 - run_SNPslice: A function created by Kathryn Murie, borrowed here
 to run the negative binomial model of SNP-Slice.
 - infer_SNPslice_freqs: A second function created by Kathryn Murie to estimate
 the frequencies of genotypes by counting the number of samples that contain
 each genotype divided by the total number of genotype observations in the
 population. Modified slightly to output the number of genotype observations in
 the population.
 - get_alternates_df: Reads the table created by create_ref_alt into an R data
 frame
 - genotype_to_stave: takes genotypes formatted as 1's and 0's, reference and
 alternate alleles, and original locus amino acid names, and reformats genotypes
 into STAVE format
 - lookup_genotype - looks up whether each locus within the genotype should be
 listed as reference or alternate, as well as the amino acid associated with
 each
 - genotypes_from_freqs - a wrapper function that creates final STAVE formatted
 output datasets
 - subset_groups - main running of the processing pipeline by group
