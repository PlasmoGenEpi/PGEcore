# Add Reference Sequences to BED Table

This folder contains two scripts to add reference sequences to a BED file containing target information. The first method uses a FASTA file where each target has the reference for that target already extracted. The second takes the full genome and extracts the reference for each target using the coordinates in the BED file. 

Contents: 
* [add\_ref\_seqs\_with\_fasta.R](#add_ref_seqs_with_fastar)
	* [Tool Information](#tool-information)
	* [Script Usage](#script-usage)
* [add\_ref\_seqs\_with\_genome.R](#add_ref_seqs_with_genomer)
	* [Tool Information](#tool-information-1)
	* [Script Usage](#script-usage-1)
## add\_ref\_seqs\_with\_fasta.R 

### Tool Information

Take a bed file containing panel location information and add reference sequences from a fasta file with the corresponding targeted reference sequences.


## Script Usage 


```bash
Usage: ./add_ref_seqs_with_targeted_ref_fasta.R [options]


Options:
	--ref_bed=REF_BED
		a bed file containing the reference location of the ref_seq, columns should be #chrom, start, end, target_name, length, strand

	--target_fasta=FASTA
		a fasta file with the ref sequences for the targets, the names of the records should match up with the target_name of the --ref_bed file

	--out=OUT
		the out file to write to

	--overwrite
		overwrite the output if it already exists

	-h, --help
		Show this help message and exit
```

An example BED file for `--ref_bed` can be found [here](../../data/example_panel_info.bed) and an example FASTA file for `--target_fasta` can be found [here](../../data/example_PMO_insert_locs_of_panel_refseqs.fasta) 

You can test the script with the example data using the following command from within this folder.

```
Rscript add_ref_seqs_with_targeted_ref_fasta.R --ref_bed ../../data/example_panel_info.bed --target_fasta ../../data/example_PMO_insert_locs_of_panel_refseqs.fasta --out example_panel_info_with_ref.bed
```

# add\_ref\_seqs\_with\_genome.R 


## Tool Information

Take a bed file containing panel location information and extract the corresponding reference sequences from a full genome fasta file, adding them to the bed file.

## Script Usage


```bash
Usage: ./add_ref_seqs_with_full_genome_ref_fasta.R [options]


Options:
	--ref_bed=REF_BED
		a bed file containing the reference location of the ref_seq, columns should be #chrom, start, end, target_name, length, strand

	--genome_fasta=GENOME
		a genome file to extract the ref_seq

	--out=OUT
		the out file to write to

	--overwrite
		overwrite the output if it already exists

	-h, --help
		Show this help message and exit
```

An example BED file for `--ref_bed` can be found [here](../../data/example_panel_info.bed). We don't include an example full genome FASTA file for `--genome_fasta` as it would be too large, but an example can be found [here](../../data/example_PMO_insert_locs_of_panel_refseqs.fasta) (TODO: come back and add link to place to download reference seq for full genome)

After you have downloaded an example full genome reference you can test the script by running a command similar to the below from within this folder.

```
Rscript add_ref_seqs_with_full_genome_ref_fasta.R --ref_bed  ../../data/example_panel_info.bed --genome_fasta Pf3D7.fasta --out example_panel_info_with_ref.bed

```




