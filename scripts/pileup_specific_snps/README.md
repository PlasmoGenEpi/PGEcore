# Pileup SNPs from microhaplotype data

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
		a bed file containing the snps of interest to calculate pileup, is genomic location of the snp position so all snps should be of size 1, should have columns, #chrom, start, end, name, length, strand

	--output_directory=OUTPUT_DIRECTORY
		the output directory to write results to

	--select_target_ids=SELECT_TARGET_IDS
		only process these targets

	--select_specimen_ids=SELECT_SPECIMEN_IDS
		only process these samples

	--overwrite_dir
		overwrite the output directory if it already exists

	--collapse_calls_by_summing
		collapse SNPs calls by summing across targets, by default the target with the highest read count will be used as the final call

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

The output will have 3-4 files. 

*  [snp\_calls.tsv.gz](#snp_callstsvgz)
*  [collapsed\_snp\_calls.tsv.gz](#collapsed_snp_callstsvgz)  
*  [snps\_covered\_by\_target\_samples\_info.tsv](snps_covered_by_target_samples_infotsv)  

#### snp\_calls.tsv.gz 

This has the snps calls per loci per seq per sample, this all calls even if the input targets overlap for SNPs so SNPs could occur more than once in the file if covered by multiple targets  

*  specimen_id - The name of the specimen/sample from the allele table  
*  target_id - The name of the target  
*  read_count - The read count for this sample for this SNP  
*  seq - the microhaplotype sequence  
*  chrom - the chromosome of the SNP
*  pos - the chromosome position of the SNP (0-based)  
*  snp_name - the name from the input file **\-\-ref\_bed**  input 
*  ref_base - base of the reference at this position 
*  seq_base - the base of the micorhaplotype sequence at this position  
*  isBiallelic - TRUE or FALSE for whether this SNP within this allele table was biallelic or not  

Example

|specimen\_id|target\_id|read\_count|seq|chrom|pos|snp\_name|ref\_base|seq\_base|isBiallelic|  
|:---|:---|:---|:---|:---|:---|:---|:---|:---|:---|  
|Laos2017-01|Pf3D7\_01\_v3-0181544-0181729|4648|TTTCATTATTGTTTTCATTCTTTTTTTAACGAAAACTATTCATCTCAAAAATATAAGATATTTTATATGACGAATGCCATTGTATTTTTTGTTACGTAAAACCTGACTTCTTCAGGGAAAACACATGCGCATTTTCACCAATTTTTGCCTAAGCTTATTATAAAAAGTATATTAAATGTATGACT|Pf3D7\_01\_v3|181574|Pf3D7\_01\_v3-181574-181575|G|G|TRUE|
|Laos2017-01|Pf3D7\_01\_v3-0181544-0181729|4648|TTTCATTATTGTTTTCATTCTTTTTTTAACGAAAACTATTCATCTCAAAAATATAAGATATTTTATATGACGAATGCCATTGTATTTTTTGTTACGTAAAACCTGACTTCTTCAGGGAAAACACATGCGCATTTTCACCAATTTTTGCCTAAGCTTATTATAAAAAGTATATTAAATGTATGACT|Pf3D7\_01\_v3|181658|Pf3D7\_01\_v3-181658-181659|A|G|TRUE|
|Laos2017-01|Pf3D7\_01\_v3-0528889-0529073|5032|ATTTGATTCTTTTTAATGAAAAAGAAGCTAAAGATATGTCAGACGATATAATTTCCCAACAAAAACGTTATTGCTCTACCAATATTCATAGTAATTATAATAATAAAATATGTATATGTAAAAATAAGCGACATCATAACAAAAGAGGGAAAGGAATAAAGCATCCTGACATACATCAAAAGGA|Pf3D7\_01\_v3|528906|Pf3D7\_01\_v3-528906-528907|G|G|TRUE|
|Laos2017-01|Pf3D7\_01\_v3-0528889-0529073|5032|ATTTGATTCTTTTTAATGAAAAAGAAGCTAAAGATATGTCAGACGATATAATTTCCCAACAAAAACGTTATTGCTCTACCAATATTCATAGTAATTATAATAATAAAATATGTATATGTAAAAATAAGCGACATCATAACAAAAGAGGGAAAGGAATAAAGCATCCTGACATACATCAAAAGGA|Pf3D7\_01\_v3|528933|Pf3D7\_01_v3-528933-528934|C|G|FALSE|

#### collapsed\_snp\_calls.tsv.gz

In some cases SNPs will be covered by multiple targets, this table will collapse them to a single call per sample. The collapsing is handled by the argument `--collapse_calls_by_summing`, the default is to pick the 'best' target which is the one with the highest read count, if `--collapse_calls_by_summing` is turned on, the read counts are instead summed accross targets. The columns in this table are similar to **snp\_calls.tsv.gz** above and will contain the columns specimen\_id, target\_id, chrom, pos, snp\_name, ref\_base, seq\_base, read\_count, isBiallelic but in addition will have these two columns 

*  best\_target\_id	- The 'best' target that had the highest read coverage for this SNP 
*  covered\_by\_target\_ids - a comma separated list of all the targets that cover this SNP  

Example

|specimen\_id|target\_id|chrom|pos|snp\_name|ref\_base|seq\_base|best\_target\_id|covered\_by\_target\_ids|read\_count|isBiallelic|
|:---|:---|:---|:---|:---|:---|:---|:---|:---|:---|:---|
|Laos2017-01|Pf3D7\_01\_v3-0181544-0181729|Pf3D7\_01\_v3|181574|Pf3D7\_01\_v3-181574-181575|G|G|Pf3D7\_01\_v3-0181544-0181729|Pf3D7\_01\_v3-0181544-0181729|4648|TRUE|
|Laos2017-01|Pf3D7\_01\_v3-0181544-0181729|Pf3D7\_01\_v3|181658|Pf3D7\_01\_v3-181658-181659|A|G|Pf3D7\_01\_v3-0181544-0181729|Pf3D7\_01\_v3-0181544-0181729|4648|TRUE|
|Laos2017-01|Pf3D7\_07\_v3-0403499-0403683|Pf3D7\_07\_v3|403611|Pf3D7\_07\_v3-403611-403612|T|T|Pf3D7\_07\_v3-0403499-0403683|Pf3D7\_07\_v3-0403499-0403683,Pf3D7\_07\_v3-0403507-0403717|5884|TRUE|
|Laos2017-01|Pf3D7\_07\_v3-0403499-0403683|Pf3D7\_07\_v3|403612|Pf3D7\_07\_v3-403612-403613|G|G|Pf3D7\_07\_v3-0403499-0403683|Pf3D7\_07\_v3-0403499-0403683,Pf3D7\_07\_v3-0403507-0403717|5884|TRUE|


#### snps\_covered\_by\_target\_samples\_info.tsv 

This has info about which SNPs you supplied are actual covered and by how many samples within the allele table 

Has the following columns, the first 6 of which are from the input `--ref_bed ` file and are bed standard columns (\#chrom, start, end, name, length, strand) and include 

*  covered_by_target - which targets within the allele table covered this SNP 
*  chrom - chromosome again 
*  pos - zero based position of the SNP 
*  ref_base - the base of the reference at this SNP  
*  n_samples - the number of samples covering this SNP  
*  total_samples - the total number of samples within the allele table 

Example: 

|#chrom|start|end|name|length|strand|covered\_by\_target|chrom|pos|ref\_base|n\_samples|total\_samples|
|:---|:---|:---|:---|:---|:---|:---|:---|:---|:---|:---|:---|
|Pf3D7\_01\_v3|145514|145515|Pf3D7\_01\_v3-145514-145515|1|+|uncovered|NA|NA|NA|NA|NA|
|Pf3D7\_01\_v3|145559|145560|Pf3D7\_01\_v3-145559-145560|1|+|uncovered|NA|NA|NA|NA|NA|
|Pf3D7\_01\_v3|181574|181575|Pf3D7\_01\_v3-181574-181575|1|+|Pf3D7\_01\_v3-0181544-0181729|Pf3D7\_01\_v3|181574|G|25|25|
|Pf3D7\_01\_v3|181658|181659|Pf3D7\_01\_v3-181658-181659|1|+|Pf3D7\_01\_v3-0181544-0181729|Pf3D7\_01\_v3|181658|A|25|25|

### Example Script Usage

Using the files in the github `data/` folder and run within the same folder of `pileup_specific_snps.R`

```
./pileup_specific_snps.R --allele_table ../../data/example2_allele_table.tsv --ref_bed ../../data/example_PMO_insert_locs_of_panel.bed --output_directory testing --snps_of_interest ../../data/MAD4HATTER_coveredSnps.bed  --overwrite_dir
```
