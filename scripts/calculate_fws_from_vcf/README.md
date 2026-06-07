# calculate\_fws\_from\_vcf.R

Contents:
* [Tool Information](#tool-information)
* [Installation](#installation)
* [Script Usage](#script-usage)


## Tool Information

Calculate the `Fws` within-host diversity statistic for each sample from a VCF
using `moimix::getFws()`. `Fws` ranges from 0 (high within-host diversity) to 1
(clonal/monoclonal); it is computed from the per-sample allelic depths
(`FORMAT/AD`), so the input VCF must carry `AD` -- e.g. one built by
`snp_calls_to_vcf.R`.

The VCF is first converted to a GDS (`SeqArray::seqVCF2GDS`). The GDS is only
(re)created when it is missing or older than the VCF (use `--force` to always
rebuild). By default the GDS path is derived from the input by swapping the
`.vcf`/`.vcf.gz` suffix for `.gds`. The output is a TSV of `specimen_name` and
`fws`, sorted by `fws`; pass `--population_name` to tag every row with a
population label for later merging across populations.


## Installation

`moimix` is not on CRAN/Bioconductor, so install it (and its dependencies) from
within R:

```r
setRepositories(ind = 1:3)
install.packages(c("devtools", "remotes", "SeqArray", "optparse", "dplyr", "readr", "tidyr"))
install.packages(c("SeqVarTools", "BiocParallel"))
remotes::install_github("bahlolab/moimix")
```


## Script Usage

```bash
Usage: ./calculate_fws_from_vcf.R -i input.vcf.gz -o fws_result.tsv [options]

Calculate Fws from a VCF using moimix::getFws().

Options:
	-i FILE, --input=FILE
		Input VCF file (.vcf or .vcf.gz) [required]

	-o FILE, --output=FILE
		Output TSV of Fws results [default: fws_result.tsv]

	-g FILE, --gds=FILE
		GDS file path. If omitted, it is derived from the input VCF by replacing the .vcf/.vcf.gz suffix with .gds

	-p NAME, --population_name=NAME
		Optional population name. If provided, adds a 'population_name' column with this value to the output (useful for combining multiple populations later)

	-f, --force
		Force re-creation of the GDS even if an up-to-date one exists

	-v, --verbose
		Print progress messages (otherwise run silently)

	-h, --help
		Show this help message and exit
```


```
./scripts/calculate_fws_from_vcf/calculate_fws_from_vcf.R -i data/example_snp_calls.vcf -o example_fws_result.tsv
```
