# add\_ref\_seqs\_with\_fasta.R 

Contents: 
* [Tool Information](#tool-information)
* [Script Usage](#script-usage)


## Tool Information

Take a bed file containing panel location information and add reference sequences from a fasta file with the corresponding targeted reference sequences.


## Script Usage 


```bash
Usage: ./add_ref_seqs_with_targeted_ref_fasta.R [options]


Options:
	--ref_bed=REF_BED
		a bed file containing the reference location of the ref_seq, no column names but the first 6 columns should be chrom, start, end, target_id, length, strand

	--target_fasta=FASTA
		a fasta file with the ref sequences for the targets, the names of the records should match up with the target_id of the --ref_bed file

	--out=OUT
		the out file to write to

	--overwrite
		overwrite the output if it already exists

	-h, --help
		Show this help message and exit
```


```
./add_ref_seqs_with_targeted_ref_fasta.R --ref_bed example_PMO_insert_locs_of_panel.bed --target_fasta example_PMO_insert_locs_of_panel_refseqs.fasta --out example_PMO_insert_locs_of_panel_withRefSeqs.bed
```


# add\_ref\_seqs\_with\_genome.R 


## Tool Information - genome

Take a bed file containing panel location information and extract the corresponding reference sequences from a full genome fasta file, adding them to the bed file.

## Script Usage - genome


```bash
Usage: ./add_ref_seqs_with_full_genome_ref_fasta.R [options]


Options:
	--ref_bed=REF_BED
		a bed file containing the reference location of the ref_seq, no column names but the first 6 columns should be chrom, start, end, target_id, length, strand

	--genome_fasta=GENOME
		a genome file to extract the ref_seq

	--out=OUT
		the out file to write to

	--overwrite
		overwrite the output if it already exists

	-h, --help
		Show this help message and exit
```


```
./add_ref_seqs_with_full_genome_ref_fasta.R --ref_bed  example_PMO_insert_locs_of_panel.bed --genome_fasta Pf3D7.fasta --out example_PMO_insert_locs_of_panel_withRefSeqs.bed

```




