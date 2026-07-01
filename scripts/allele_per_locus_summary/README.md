# Allele per Locus Summary

Contents: 
* [Tool Information](#tool-information)
* [Script Usage](#script-usage)

## Tool Information

This tool provides a simple count of alleles in the dataset. By locus, it counts:
* Total Allele Count: The total number of alleles (length(seq)).
* Unique Allele Count: The count of unique alleles (length(unique(seq))).
* Allele Singlets: The number of alleles that appear only once sum(table(seq) == 1)  

## Script Usage 

```bash
Usage: allele_per_locus_summary.R [options]


Options:
        --allele_table=ALLELE_TABLE
                TSV containing allele present/absent per specimen, with the
       columns: specimen_name, target_name, seq

        -h, --help
                Show this help message and exit
```

An example allele table can be found [here](../../data/example_allele_table.tsv). In this example, the seq column contains the sequence of the allele. For this script to work, that column can include any unique identifier for the allele. Identifiers need to be consistent across samples (i.e. the ID for the same allele must be the same in all samples.)

You can use the below command to test running this script on example data within this repository 
```
Rscript allele_per_locus_summary.R --allele_table ../../data/example_allele_table.tsv
```

TODO: this script shuld allow the user to change the name of the output file.