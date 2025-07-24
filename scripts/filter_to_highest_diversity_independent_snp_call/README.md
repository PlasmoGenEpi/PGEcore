# Filter SNP calls to relatively independent SNPs

Contents: 

* [Tool Information](#tool-information)
* [Script Usage](#script-usage)


## Tool Information
This will take a SNP pileup (performed by `pileup_specific_snps.R`) and filter to the highest diveristy (by expected heterozygosity(he)) within a specific distance between SNPs (default: 10,000bp and can be changed with `--mindist_between_snps`). 


## Script Usage 


```bash
Usage: ./filter_to_highest_diversity_independent_snp_call.R [options]


Options:
	--snp_table=SNP_TABLE
		TSV containing at least the columns: specimen_id, target_id, chrom, pos, snp_name, ref_base, seq_base, read_count, is_biallelic

	--output_fnp=OUTPUT_FNP
		the output filename path for the output filtered file

	--mindist_between_snps=MINDIST_BETWEEN_SNPS
		the minimum distance between SNPs to consider them 'independent'. Default: 10000

	--select_target_ids=SELECT_TARGET_IDS
		only process these targets

	--select_specimen_ids=SELECT_SPECIMEN_IDS
		only process these samples

	--overwrite
		overwrite the output fnp if it already exists. Default: FALSE

	--only_biallelic
		filter to only biallelic SNPs too. Default: FALSE

	--only_informative
		filter off SNPs with no variation. Default: FALSE

	-h, --help
		Show this help message and exit
```

Two required options of `--snp_table` and `--output_fnp `. Details below for `--snp_table`. 


### \-\-snp_table 
The expected input is the collapsed SNP calls from `pileup_specific_snps/pileup_specific_snps.R` though any table with columns of `specimen_id, target_id, chrom, pos, snp_name, ref_base, seq_base, read_count, is_biallelic` will work as long as the input has only 1 call per `specimen_id` and `seq_base` (e.g. SNP calls that have the same SNPs accross different target will cause the script to throw an error). 

#### Example

Here is example SNP table that could be input into the script. 

|specimen\_id|target\_id|chrom|pos|snp\_name|ref\_base|seq\_base|best\_target\_id|covered\_by\_target\_ids|read\_count|isBiallelic|
|:---|:---|:---|:---|:---|:---|:---|:---|:---|:---|:---|
|Laos2017-01|Pf3D7\_01\_v3-0181544-0181729|Pf3D7\_01\_v3|181574|Pf3D7\_01\_v3-181574-181575|G|G|Pf3D7\_01\_v3-0181544-0181729|Pf3D7\_01\_v3-0181544-0181729|4648|TRUE|
|Laos2017-01|Pf3D7\_01\_v3-0181544-0181729|Pf3D7\_01\_v3|181658|Pf3D7\_01\_v3-181658-181659|A|G|Pf3D7\_01\_v3-0181544-0181729|Pf3D7\_01\_v3-0181544-0181729|4648|TRUE|
|Laos2017-01|Pf3D7\_07\_v3-0403499-0403683|Pf3D7\_07\_v3|403611|Pf3D7\_07\_v3-403611-403612|T|T|Pf3D7\_07\_v3-0403499-0403683|Pf3D7\_07\_v3-0403499-0403683,Pf3D7\_07\_v3-0403507-0403717|5884|TRUE|
|Laos2017-01|Pf3D7\_07\_v3-0403499-0403683|Pf3D7\_07\_v3|403612|Pf3D7\_07\_v3-403612-403613|G|G|Pf3D7\_07\_v3-0403499-0403683|Pf3D7\_07\_v3-0403499-0403683,Pf3D7\_07\_v3-0403507-0403717|5884|TRUE|

### Output

The output will be a filtered table of the input so it will look the same as the input table. 

### Example Script Usage

Using the files in the github `data/` folder and run within the same folder of `filter_to_highest_diversity_independent_snp_call.R`

```bash
# run a pileup of SNPs 
./../pileup_specific_snps/pileup_specific_snps.R --allele_table ../../data/example2_allele_table.tsv --ref_bed ../../data/example_PMO_insert_locs_of_panel.bed --output_directory testing --snps_of_interest ../../data/MAD4HATTER_coveredSnps.bed  --overwrite_dir

# filter SNPS 
./filter_to_highest_diversity_independent_snp_call.R --snp_table testing/collapsed_snp_calls.tsv.gz --out independent_collapsed_snp_calls.tsv.gz

## use only biallelic and only informative snps 
./filter_to_highest_diversity_independent_snp_call.R --snp_table testing/collapsed_snp_calls.tsv.gz --out independent_collapsed_snp_calls.tsv.gz --only_biallelic --only_informative --overwrite

```

