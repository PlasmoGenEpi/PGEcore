# Insert tool/ code title

Contents: 
* [Tool Information](#tool-information)
* [Script Usage](#script-usage)


## Tool Information
This will perform a pileup of specific snp positions from the microhaplotypes


## Script Usage 


```bash
Usage: ./pileup_specific_snps.R [options]


Options:
	--allele_table=ALLELE_TABLE
		TSV containing the columns: specimen_id, target_id, read_count, seq

	--ref_bed=REF_BED
		a bed file containing the reference location of the ref_seq, no column names but the first 6 columns should be chrom, start, end, target_id, length, strand, ref_seq

	--snps_of_interest=SNPS_OF_INTEREST
		a bed file containing the snps of interest to calculate pileup, is genomic location of the snp position so all loci should be of size 1, should have columns, #chrom, start, end, name, length, strand

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


Four required options: `--allele_table`, `--ref_bed`, `--snps_of_interest`, `--output_directory`. Options are available to handling which targets and/or specimens to analyze. 

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

### --snps\_of\_interest

Supplying 1 base loci of the snps positions, requires #chrom, start, end, name, length, strand columns. The gene column can be free form and gene_id should be a canonical identifier. 

|#chrom|start|end|name|length|strand|
|---|---|---|---|---|---|
|Pf3D7\_01\_v3|145514|145515|Pf3D7\_01\_v3-145514-145515|1|+|
|Pf3D7\_01\_v3|145559|145560|Pf3D7\_01\_v3-145559-145560|1|+|
|Pf3D7\_01\_v3|162980|162981|Pf3D7\_01\_v3-162980-162981|1|+|
|Pf3D7\_01\_v3|163070|163071|Pf3D7\_01\_v3-163070-163071|1|+|
|Pf3D7\_01\_v3|163079|163080|Pf3D7\_01\_v3-163079-163080|1|+|



### --output\_directory

An output directory, can be overwritten with `--overwrite_dir`



### Example 

Using the files in the github `data/` folder

```
./pileup_specific_snps.R --allele_table ../../data/example2_allele_table.tsv --ref_bed ../../data/example_PMO_insert_locs_of_panel.bed --output_directory testing --snps_of_interest ../../data/MAD4HATTER_coveredSnps.bed  --overwrite_dir
```
