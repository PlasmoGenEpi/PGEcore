# add\_ref\_seqs\_with\_fasta.R 

Contents: 
* [Tool Information](#tool-information)
* [Script Usage](#script-usage)


## Tool Information

Take a bed file with location of panel info and add reference sequences from a fasta file to that file. 



## Script Usage 


```bash
Usage: ./add_ref_seqs_with_fasta.R [options]


Options:
	--ref_bed=REF_BED
		a bed file containing the reference location of the ref_seq, no column names but the first 6 columns should be chrom, start, end, target_id, length, strand

	--fasta=FASTA
		a fasta file with the ref sequences, the names of the records should match up with the target_id of the --ref_bed file

	--out=OUT
		the out file to write to

	--overwrite
		overwrite the output if it already exists

	-h, --help
		Show this help message and exit
```


```
./add_ref_seqs_with_fasta.R --ref_bed example_PMO_insert_locs_of_panel.bed --fasta example_PMO_insert_locs_of_panel_refseqs.fasta --out example_PMO_insert_locs_of_panel_withRefSeqs.bed
```


# add\_ref\_seqs\_with\_genome.R 


## Tool Information - genome

Take a bed file with location of panel info and add reference sequences from a genome file to that file. 


## Script Usage - genome


```bash
Usage: ./add_ref_seqs_with_genome.R [options]


Options:
	--ref_bed=REF_BED
		a bed file containing the reference location of the ref_seq, no column names but the first 6 columns should be chrom, start, end, target_id, length, strand

	--genome=GENOME
		a genome file to extract the ref_seq

	--out=OUT
		the out file to write to

	--overwrite
		overwrite the output if it already exists

	-h, --help
		Show this help message and exit
```


```
./add_ref_seqs_with_genome.R --ref_bed  example_PMO_insert_locs_of_panel.bed --genome Pf3D7.fasta --out example_PMO_insert_locs_of_panel_withRefSeqs.bed

```




