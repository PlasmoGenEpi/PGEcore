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
To run on all targets (write files to current working directory):
```
Rscript scripts/malariaem_wrapper/malariaem_wrapper.R --allele_table \
    data/example2_allele_table.tsv \
    --subset_targets FALSE \
    --coi_range "1, 2, 3, 4" \
```

To run on all targets (write files to specific directory, *note: this must be created beforehand if it doesn't exist already*):
```
Rscript scripts/malariaem_wrapper/malariaem_wrapper.R --allele_table \
    data/example2_allele_table.tsv \
    --subset_targets FALSE \
    --coi_range "1, 2, 3, 4" \
    --output_directory results_malariaem/
```

To run with subsetting by specific groups (recommended if COI and number of targets is large since `malaria.em` is slow):
```
Rscript scripts/malariaem_wrapper/malariaem_wrapper.R \
  --allele_table data/example2_allele_table.tsv \
  --subset_targets TRUE \
  --target_groups data/example_target_groups.tsv \
  --coi_range "1, 2, 3, 4" \
```