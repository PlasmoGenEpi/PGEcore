# SNP-Slice Wrapper

Contents: 
* [Tool Information](#tool-information)
* [Script Usage](#script-usage)

## Tool Information

SNP-Slice is a Bayesian nonparametric method for resolving multi-strain 
infections using slice sampling. The algorithm simultaneously estimates strain 
haplotypes and links them to hosts from SNP data.

This wrapper runs SNP-Slice and calculates population-level multilocus allele 
frequencies (MLAF) and sample-level complexity of infection (COI) estimates from 
the SNP-Slice results. This functionality is supported by the 
[snp.slicer](https://github.com/PlasmoGenEpi/snp.slicer) R package, which is a 
repackage of the original SNP-Slice code.

### Existing Resources

SNP-Slice is described in: Ju, N., Liu, J., & He, Q. (2024). SNP-slice resolves 
mixed infections: Simultaneously unveiling strain haplotypes and linking them to 
hosts. Bioinformatics (Oxford, England), 40(6), btae344. 
https://doi.org/10.1093/bioinformatics/btae344

The original SNP-Slice code can be found 
[here](https://github.com/nianqiaoju/snp-slice) and the snp.slicer repackaging 
of this code is [here](https://github.com/PlasmoGenEpi/snp.slicer).

## Inputs and Outputs

The script requires an allele table (`--allele_table`) and a loci groups table 
(`--loci_groups_input`). The allele table should only contain biallelic 
genotypes. A variety of optional inputs specify the column names in these tables 
and parameters for running SNP-Slice. For specifics on input formats, refer to 
the argument documentation in the script.

The MLAF output will be a TSV with group\_id, variant, and freq columns. If the 
`--use_mcmc_for_af_and_coi` flag was not used, it will also have allele\_count 
and allele\_total columns.

The COI output will be a TSV with a column matching `--specimen_name_col` and a 
coi column. If the `--use_mcmc_for_af_and_coi` flag was used, it will also have 
coi\_sd, coi\_lower, and coi\_upper columns. The latter three are the standard 
deviation and lower and upper bounds of the 95% credible interval of the COI 
estimate.

## Script Usage

```{r}
# Basic usage
scripts/snpslice_wrapper/snpslice_wrapper.R \
    --allele_table data/example_amino_acid_calls.tsv \
    --loci_groups_input data/example_loci_groups.tsv \
    --mlaf_output mlaf.tsv \
    --coi_output coi.tsv

# Use MCMC results for estimating MLAF and COI
scripts/snpslice_wrapper/snpslice_wrapper.R \
    --allele_table data/example_amino_acid_calls.tsv \
    --loci_groups_input data/example_loci_groups.tsv \
    --mlaf_output mlaf.tsv \
    --coi_output coi.tsv \
    --use_mcmc_for_af_and_coi
```
