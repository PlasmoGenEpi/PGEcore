# Translating amino acid loci of interest

Contents: 
* [Tool Information](#tool-information)
* [Script Usage](#script-usage)




## Tool Information

This script handles translating the microhaplotype sequences into loci of interest, commonly known drug resistance loci. 

This takes in an allele table, a reference sequence table to align to, the bed location of those alleles, and a bed of the loci of interest. 


## Script Usage 

requires the below packages 

```r
BiocManager::install("pwalign")
BiocManager::install("Biostrings")
install.packages("tibble", "dplyr", "stringr", "readr", "optparse")
```


```
Options:
	--allele_table=ALLELE_TABLE
		TSV containing the columns: specimen_id, target_id, read_count, seq

	--ref_bed=REF_BED
		a bed file containing the reference location of the ref_seq, no column names but the first 6 columns should be chrom, start, end, target\_id, length, strand, ref\_seq

	--loci_of_interest=LOCI_OF_INTEREST
		a bed file containing the loci of interest to translate, is genomic location of the codon position so all loci should be of size 3, should correspond to the same, should have columns, #chrom, start, end, name, length, strand, gene, gene_id, aa_position 

	--output_directory=OUTPUT_DIRECTORY
		the output directory to write results to

	--select_target_ids=SELECT_TARGET_IDS
		only process these targets

	--select_specimen_ids=SELECT_SPECIMEN_IDS
		only process these samples

	--overwrite_dir
		overwrite the output directory if it already exists

	--collapse_calls_by_summing
		collapse amino acid calls by summing across targets, by default the target with the highest read count will be used as the final call

	-h, --help
		Show this help message and exit
```

Four required options: `--allele_table`, `--ref_bed`, `--loci_of_interest`, `--output_directory`. Options are available to handling which targets and/or specimens to analyze. 

### --allele_table

Table of microhaplotype calls, needs 4 columns named: 

| specimen_id | target_id | read_count | seq |  
| --- | --- | --- | --- |  
| Laos2017-01 | Pf3D7\_01\_v3-0181544-0181729 | 4648.0 |   TTTCATTATTGTTTTCATTCTTTTTTTAACGAAAACTATTCATCTCAAAAATATAAGATATTTTATATGACGAATGCCATTGTATTTTTTGTTACGTAAAACCTGACTTCTTCAGGGAAAACACATGCGCATTTTCACCAATTTTTGCCTAAGCTTATTATAAAAAGTATATTAAATGTATGACT |  
|Laos2017-01|Pf3D7\_01\_v3-0528889-0529073|5032.0|  ATTTGATTCTTTTTAATGAAAAAGAAGCTAAAGATATGTCAGACGATATAATTTCCCAACAAAAACGTTATTGCTCTACCAATATTCATAGTAATTATAATAATAAAATATGTATATGTAAAAATAAGCGACATCATAACAAAAGAGGGAAAGGAATAAAGCATCCTGACATACATCAAAAGGA|
|Laos2017-01|Pf3D7\_02\_v3-0290610-0290789|5220.0|  TTTATTATAAGGATTATAATTATTATCATTCTTATTATATGGATTAATATAATTTCGATCATAATCATTACTGTTTTTATCATATTGTGGATAATTATTTCTATTTCTATTTTGATTAAAATCCGATTGTCTATATCTAACTCTAATATCTGTAGATGCTAAATCATTAGGTAATGGGA|  


### --ref\_bed

This supplies the location of the panel and the reference sequence of that panel 

Needs columns "#chrom", "start", "end", "target\_id", "length", "strand", "ref\_seq"

|#chrom|start|end|target\_id|length|strand|ref\_seq|
|---|---|---|---|---|---|---|
|Pf3D7\_01\_v3|181544|181729|Pf3D7\_01\_v3-0181544-0181729|185|+|TTTCATTATTGTTTTCATTCTTTTTTTAACGAAAACTATTCATCTCAAAAATATAAGATATTTTATATGACGAATGCCATTGTATTTTTTGTTACGTAAAACCTGACTTCTTCAAGGAAAACACATGCGCATTTTCACCAATTTTTGCCTAAGCTTATTATAAAAAGTATATTAAATGTATGACT|
|Pf3D7\_01\_v3|528889|529073|Pf3D7\_01\_v3-0528889-0529073|184|+|ATTTGATTCTTTTTAATGAAAAAGAAGCTAAAGATATGTCAGACCATATAATTTCCCAACAAAAACGTTATTGCTCTACCAATATTCATAGTAATTATAATAATAAAATATGTATATGTAAAAATAAGCGACATCATAACAAAAGAGGGAAAGGAATAAAGCATCCTGACATACATCAAAAGGA|
|Pf3D7\_02\_v3|290610|290789|Pf3D7\_02\_v3-0290610-0290789|179|+|TTTATTATAAGGATTATAATTATTATCATTCTTATTATATGGATTAATATAATTTCGATCATTATCATTACTGTTTTTATCATATTGTGGATAATTATTTCTATTTCTATTTTGATTAAAATCCGATTGTCTATATCTAACTCTAATATCTGTAGATGCTAAATCATTAGGTAATGGGA|

### --loci\_of\_interest

Supplying 3 base loci of the codon positions of the amino acid positions, requires #chrom, start, end, name, length, strand, gene, gene\_id, aa\_position columns. The gene column can be free form and gene_id should be a canonical identifier. 

|#chrom|start|end|name|length|strand|gene|aa\_position|gene\_id|
|---|---|---|---|---|---|---|---|---|
|Pf3D7\_01\_v3|268384|268387|PF3D7\_0106300.1-AA263|3|-|atp6|263|PF3D7\_0106300|
|Pf3D7\_01\_v3|267880|267883|PF3D7\_0106300.1-AA431|3|-|atp6|431|PF3D7\_0106300|
|Pf3D7\_01\_v3|267304|267307|PF3D7\_0106300.1-AA623|3|-|atp6|623|PF3D7\_0106300|

### --output\_directory

An output directory, can be overwritten with `--overwrite_dir`



### Example 

Using the files in the github `data/` folder

```
../scripts/translate_loci_of_interest/translate_loci_of_interest.R --output_directory example_output --allele_table example2_allele_table.tsv  --ref_bed example_PMO_insert_locs_of_panel.bed  --loci_of_interest example_principal_resistance_marker_info_table.bed --overwrite_dir
```

