# Freq Estimation Model

Contents: 
* [Tool Information](#tool-information)
* [Script Usage](#script-usage)

## Tool Information - IF NOT ON PGEforge
If no tutorial on PGEforge then fill in the table below and delete the sections for [on pgeforge](#tool-information---if-on-pgeforge) and [bespoke code](#tool-information---if-bespoke-code).

| Tool Summary    |  |
| -------- | ------- |
| Main use-cases | Estimate multi-locus allele frequence |
| Authors | Aimee Taylor |
| Latest version | INSERT INFO |
| License | INSERT INFO |
| Website | INSERT INFO |
| Code repository | INSERT INFO |
| Publication | INSERT INFO |

### Purpose

Freq Estimation Model (FEM) estimates the population-level frequencies of multi-allele genotypes. It was developed for analyzing anti-malarial drug resistance. For the full documentation by Aimee Taylor, see the [full documentation](https://github.com/aimeertaylor/FreqEstimationModel/blob/master/README.md).

### Existing resources

- [Install instructions] https://github.com/aimeertaylor/FreqEstimationModel/blob/master/Instructions
- [Initial FEM publication] https://malariajournal.biomedcentral.com/articles/10.1186/1475-2875-13-102
- [In-depth description (Aimee's thesis)] https://github.com/aimeertaylor/FreqEstimationModel/blob/master/inst/Thesis_methods_chapter.pdf

## Script Usage 
--aa_calls "example_amino_acid_calls.tsv" --coi "example_coi_table.tsv" --groups "example_loci_groups.tsv" --mlaf "output"