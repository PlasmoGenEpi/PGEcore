# Translating amino acid loci of interest

Contents: 
* [Tool Information](#tool-information)
* [Script Usage](#script-usage)




## Tool Information

This script handles translating the microhaplotype sequences into loci of interest, commonly known drug resistance loci. 

This takes in an allele table, a reference sequence table to align to, the bed location of those alleles, and a bed of the loci of interest. 


## Script Usage 

```
Options:
	--allele_table=ALLELE_TABLE
		TSV containing the columns: specimen_id, target_id, read_count, seq

	--ref_seq_table=REF_SEQ_TABLE
		TSV containing the columns: target_id, ref_seq

	--ref_bed=REF_BED
		a bed file containing the reference location of the ref_seq, no column names but the first 6 columns should be chrom, start, end, name, length, strand

	--loci_of_interest=LOCI_OF_INTEREST
		a bed file containing the loci of interest to translate, is genomic location of the codon position so all loci should be of size 3, should correspond to the same 

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

Five required options: `--allele_table`, `--ref_seq_table`, `--ref_bed`, `--loci_of_interest`, `--output_directory`. Options are available to handling which targets and/or specimens to analyze. 

### --allele_table

Table of microhaplotype calls, needs 4 columns named: 

| specimen_id | target_id | read_count | seq |  
| --- | --- | --- | --- |  
| Laos2017-01 | Pf3D7\_01\_v3-0181544-0181729 | 4648.0 |   TTTCATTATTGTTTTCATTCTTTTTTTAACGAAAACTATTCATCTCAAAAATATAAGATATTTTATATGACGAATGCCATTGTATTTTTTGTTACGTAAAACCTGACTTCTTCAGGGAAAACACATGCGCATTTTCACCAATTTTTGCCTAAGCTTATTATAAAAAGTATATTAAATGTATGACT |  
|Laos2017-01|Pf3D7\_01\_v3-0528889-0529073|5032.0|  ATTTGATTCTTTTTAATGAAAAAGAAGCTAAAGATATGTCAGACGATATAATTTCCCAACAAAAACGTTATTGCTCTACCAATATTCATAGTAATTATAATAATAAAATATGTATATGTAAAAATAAGCGACATCATAACAAAAGAGGGAAAGGAATAAAGCATCCTGACATACATCAAAAGGA|
|Laos2017-01|Pf3D7\_02\_v3-0290610-0290789|5220.0|  TTTATTATAAGGATTATAATTATTATCATTCTTATTATATGGATTAATATAATTTCGATCATAATCATTACTGTTTTTATCATATTGTGGATAATTATTTCTATTTCTATTTTGATTAAAATCCGATTGTCTATATCTAACTCTAATATCTGTAGATGCTAAATCATTAGGTAATGGGA|  


### --ref\_seq\_table

Needs two columns 