# THEREALMcCOIL

Contents:

- [Tool Information](#tool-information)
- [Script Usage](#script-usage)

## Tool Information

A tutorial and information about this tool can be found [here](https://mrc-ide.github.io/PGEforge/tutorials/THEREALMcCOIL/RMCL_background.html).

### Purpose

To estimate complexity of infection and population allele frequencies used a Markov chain Monte Carlo-based method.

See the example usage for either the [categorical model](https://github.com/EPPIcenter/THEREALMcCOIL/blob/master/categorical_method/test_R_code.R) or the [proportional model](https://github.com/EPPIcenter/THEREALMcCOIL/blob/master/proportional_method/test_R_code.R).

THE REAL McCOIL assumed that different loci are independent, that different samples are independent and polygenomic infections are obtained from multiple independent infections, and that the samples were collected from a single homogeneous population.

### Existing resources

- More details about this method are available in the original [paper](https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1005348).

## Tool Information

This wrapper script simplifies the use of the THEREALMcCOIL tool by bundling all required C and R code, along with any necessary precomputed data, into a single R script. It automatically compiles the C components at runtime, reducing setup complexity and streamlining the workflow for users.

For alternative wrapper tools, see [McCOILR](https://github.com/OJWatson/McCOILR).

## Script Usage

1. Run the categorical model

```sh
Rscript ./THEREALMcCOIL_wrapper.R \
    --model categorical \
    --snp_calls_input test_data/independent_collapsed_snp_calls.tsv \
    --slaf_output test_data_output.slaf \
    --coi_output test_data_output.coi
```

2. Run the proportional model

```sh
Rscript ./THEREALMcCOIL_wrapper.R \
    --model proportional \
    --snp_calls_input test_data/independent_collapsed_snp_calls.tsv \
    --slaf_output test_data_output.slaf \
    --coi_output test_data_output.coi
```

3. Print help message.

```sh
Rscript ./THEREALMcCOIL_wrapper.R --help
```

The wrapper script provides access to all parameters required to run THEREALMcCOIL. For detailed information on available options and their usage, refer to the help message by running the script with the `--help` flag.

## Input SNP Data

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