# filter\_bialleleic\_calls.R 

Contents: 
* [Tool Information](#tool-information)
* [Script Usage](#script-usage)


## Tool Information

This will filter a table of amino acid calls and filter it to only allele calls that are bialleleic within the input population. Can also optionally write out the filtered off nonbialleleic calls 


## Script Usage 


```bash
Usage: ./filter_bialleleic_calls.R [options]


Options:
	--amino_acid_calls=AMINO_ACID_CALLS
		table with amino acid calls

	--out=OUT
		the out file to write to

	--out_nonbiallelic=OUT_NONBIALLELIC
		only process these targets

	--overwrite
		overwrite the output if it already exists

	-h, --help
		Show this help message and exit
```


### examples 


```bash
./filter_bialleleic_calls.R --amino_acid_calls ../../data/example2_amino_acid_calls.tsv --out bialleleic_example2_amino_acid_calls.tsv.tsv.gz --overwrite

./filter_bialleleic_calls.R --amino_acid_calls ../../data/example2_amino_acid_calls.tsv --out bialleleic_example2_amino_acid_calls.tsv.tsv.gz --overwrite --out_nonbiallelic nonbialleleic_example2_amino_acid_calls.tsv.tsv.gz
```

