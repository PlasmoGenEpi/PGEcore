# malaria.em

Contents: 

* [Tool Information](#tool-information)
* [Script Usage](#script-usage)

## Tool Information

### Purpose

This tool performs multi-locus haplotype frequency estimation and phasing using maximum-likelihood estimation.

### Existing resources

- The algorithm and methodology is published in [Li et al, 2007](https://doi.org/10.2202/1544-6115.1321). 

- The code base for this R package was archived on CRAN in 2012, but the code base with minor bug fixes can now be found at: https://github.com/PlasmoGenEpi/malaria.em. It can be installed via the PlasmoGenEpi R-universe at: https://plasmogenepi.r-universe.dev/malaria.em.

- Installation instructions and tutorials are on [PGEforge](https://mrc-ide.github.io/PGEforge/tutorials/malariaem/malariaem_background.html).


## Script Usage 
To run:
```
Rscript scripts/malariaem_wrapper/malariaem_wrapper.R --allele_table \
    data/example2_allele_table.tsv \
    --subset_alleles FALSE \
    --allele_names NULL \
    --coi "1, 2, 3, 4" \
    --write_results TRUE
```

To run with subsetting by specific alleles (recommended if COI and number of alleles is large since `malaria.em` is slow):
```
Rscript scripts/malariaem_wrapper/malariaem_wrapper.R \
  --allele_table data/example2_allele_table.tsv \
  --subset_alleles TRUE \
  --allele_names "Pf3D7_12_v3-0659902-0660096, Pf3D7_13_v3-1419395-1419601" \
  --coi "1, 2, 3, 4" \
  --write_results TRUE
```