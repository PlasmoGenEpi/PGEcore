# Freq Estimation Model

Contents: 
* [Tool Information](#tool-information)
* [Script Usage](#script-usage)

## Tool Information 

| Tool Summary    |  |
| -------- | ------- |
| Main use-cases | Estimate multi-locus allele frequency |
| Authors | Aimee Taylor |
| Latest version | NA |
| License | MIT |
| Website | https://github.com/aimeertaylor/FreqEstimationModel/tree/master |
| Code repository | https://github.com/aimeertaylor/FreqEstimationModel/tree/master |
| Publication | https://malariajournal.biomedcentral.com/articles/10.1186/1475-2875-13-102 |

### Purpose

FreqEstimationModel (FEM) estimates the population-level frequencies of 
multi-locus genotypes. It was developed for analyzing anti-malarial drug 
resistance. For the full documentation by  Aimee Taylor, see 
https://github.com/aimeertaylor/FreqEstimationModel/blob/master/README.md.

### Existing resources

- Install instructions: 
  https://github.com/aimeertaylor/FreqEstimationModel/blob/master/Instructions
- Initial FEM publication: 
  https://malariajournal.biomedcentral.com/articles/10.1186/1475-2875-13-102
- In-depth description (Aimee's thesis): 
  https://github.com/aimeertaylor/FreqEstimationModel/blob/master/inst/Thesis_methods_chapter.pdf

## Script Usage 

```
Rscript scripts/FreqEstimationModel_wrapper/FreqEstimationModel_wrapper.R --aa_calls data/example_amino_acid_calls.tsv --coi data/example_coi_table.tsv --groups data/example_loci_groups.tsv --mlaf_output output.tsv
```
