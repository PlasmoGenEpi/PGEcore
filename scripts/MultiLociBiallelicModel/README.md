# MultiLociBiallelicModel

Contents: 
* [Tool Information](#tool-information)
* [Script Usage](#script-usage)

## Tool Information

A tutorial and information about this tool can be found [here](https://mrc-ide.github.io/PGEforge/tutorials/MultiLociBiallelicModel/MultiLociBiallelicModel_background.html).

### Purpose

The code supplied by this paper does a maximum-likelihood (MLE) method to estimate:

haplotype frequencies and prevalence
multiplicity of infection (MOI/COI) from SNP data.
The functions here provide functionality to predict possible haplotype prevalence within the population that lead to the current set of data. Also estimates MOI/COI based on these estimates. Has to take only biallelic SNPs and compuationally can be limited by the number of loci supplied (in publication used 10 loci).

### Existing resources

- Example file can be found [here](https://github.com/Maths-against-Malaria/MultiLociBiallelicModel/blob/main/src/SNP_MLE.R)

## Script Usage 

The script can be run using the following command:

```bash
Rscript MultiLociBiallelicModel.R --aa_calls <aa_calls.tsv>
```
