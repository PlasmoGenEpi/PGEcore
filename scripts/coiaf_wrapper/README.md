# COIAF Wrapper Script

This script provides a command-line interface for estimating the Complexity of Infection (COI) using the COIAF R package. It processes SNP data and optionally population-level minor allele frequencies to estimate COI using both frequency and variant-based methods.

## Overview

The COIAF (Complexity of Infection Analysis Framework) wrapper script allows users to:
- Process SNP data from sequencing experiments
- Calculate or use provided population-level minor allele frequencies
- Estimate COI using both frequency and variant-based methods
- Output results in a standardized TSV format

## Installation

### Prerequisites

1. **R** (version 4.0 or higher)
2. **Required R packages**:
   - `coiaf` - Main COI estimation package
   - `dplyr` - Data manipulation
   - `tibble` - Data frame utilities
   - `readr` - Fast file reading/writing
   - `optparse` - Command-line argument parsing

### Installing R Packages

```r
# Install required packages
install.packages(c("dplyr", "tibble", "readr", "optparse"))

# Install COIAF (if not already installed)
install.packages('coiaf', repos = c('https://plasmogenepi.r-universe.dev', 'https://cloud.r-project.org'))
```

## Usage

### Basic Usage

```bash
Rscript coiaf_wrapper.R --snp_data <path> --output <path> [options]
```

### Required Arguments

- `--snp_data`: Path to SNP data file (TSV format)
- `--output`: Path for output file (TSV format)

### Optional Arguments

- `--plmaf`: Path to population-level minor allele frequency file (TSV format)
- `--seq_error`: Sequencing error rate (default: 0.01)
- `--max_coi`: Maximum COI to consider (default: 25)
- `--help`: Show help message

### Examples

#### Basic analysis with auto-calculated PLMAF:
```bash
Rscript coiaf_wrapper.R \
  --snp_data data/example_collapsed_snp_calls.tsv \
  --output results/coi_estimates.tsv
```

#### Analysis with custom PLMAF and parameters:
```bash
Rscript coiaf_wrapper.R \
  --snp_data data/example_collapsed_snp_calls.tsv \
  --plmaf data/example_coiaf_plmaf.tsv \
  --output results/coi_estimates.tsv \
  --seq_error 0.005 \
  --max_coi 30
```

#### Show help:
```bash
Rscript coiaf_wrapper.R --help
```

## Input File Formats

### SNP Data File (Required)

A tab-separated values (TSV) file with the following columns:

| Column | Type | Description |
|--------|------|-------------|
| `specimen_id` | String | Unique identifier for each specimen |
| `snp_name` | String | Unique identifier for each genomic target/locus |
| `read_count` | Integer | Number of reads supporting this allele |
| `seq_base` | String | Nucleotide base (A, C, G, T) |

**Example SNP data:**
```tsv
specimen_id	snp_name	read_count	seq_base
sample_001	locus_001	150	A
sample_001	locus_001	50	T
sample_001	locus_002	200	C
sample_002	locus_001	180	A
sample_002	locus_001	20	T
```

### PLMAF Data File (Optional)

A tab-separated values (TSV) file with the following columns:

| Column | Type | Description |
|--------|------|-------------|
| `snp_name` | String | Unique identifier for each genomic target/locus |
| `seq_base` | String | Nucleotide base (A, C, G, T) |
| `plmaf` | Numeric | Population-level minor allele frequency (0-1) |

**Example PLMAF data:**
```tsv
snp_name	seq_base	plmaf
locus_001	T	0.25
locus_002	G	0.15
locus_003	A	0.30
```

**Note:** If PLMAF data is not provided, the script will automatically calculate population-level minor allele frequencies from the SNP data.

## Output Format

The script outputs a tab-separated values (TSV) file with the following columns:

| Column | Type | Description |
|--------|------|-------------|
| `specimen_id` | String | Unique identifier for each specimen |
| `coi_freq` | Numeric | COI estimate using frequency method |
| `coi_variant` | Numeric | COI estimate using variant method |

**Example output:**
```tsv
specimen_id	coi_freq	coi_variant
sample_001	2.5	2.3
sample_002	1.8	1.7
sample_003	3.2	3.1
```

## Algorithm Details

The script performs the following steps:

1. **Data Processing**: Calculates within-sample minor allele frequencies (WSMAF) from read counts
2. **PLMAF Calculation**: If not provided, calculates population-level minor allele frequencies from aggregated SNP data
3. **Data Filtering**: Filters data to include only targets with valid PLMAF values
4. **COI Estimation**: Uses COIAF's `optimize_coi` function with two methods:
   - **Frequency method**: Based on allele frequency distributions
   - **Variant method**: Based on variant presence/absence patterns
5. **Output Generation**: Writes results to the specified output file

## Performance Considerations

- Processing time scales with the number of specimens and targets
- Memory usage depends on the size of input data
- For large datasets, consider processing in batches

## Troubleshooting

### Common Issues

1. **"Missing required columns" error**: Ensure your input files have the correct column names and format
2. **"File not found" error**: Check that file paths are correct and files exist
3. **"seq_error must be between 0 and 1" error**: Ensure sequencing error rate is a valid probability
4. **"max_coi must be at least 1" error**: Ensure maximum COI is a positive integer

## Citation

If you use this script in your research, please cite the COIAF package:

```
 @article{
    Paschalidis_Watson_Aydemir_Verity_Bailey_2023, 
    title={coiaf: Directly estimating complexity of infection with allele frequencies}, 
    volume={19}, 
    ISSN={1553-7358}, 
    DOI={10.1371/journal.pcbi.1010247}, 
    number={6}, 
    journal={PLOS Computational Biology}, 
    publisher={Public Library of Science}, 
    author={Paschalidis, Aris and Watson, Oliver J. and Aydemir, Ozkan and Verity, Robert and Bailey, Jeffrey A.}, 
    year={2023}, 
    month=june, 
    pages={e1010247}, 
    language={en} 
}
```

## License

This script is provided as-is for research purposes. Please refer to the COIAF package license for the underlying algorithm.
