# filter\_biallelic\_calls.R 

Contents: 
* [Tool Information](#tool-information)
* [Script Usage](#script-usage)


## Tool Information

This will filter a table of amino acid calls and filter it to only allele calls that are biallelic within the input population. Can also optionally write out the filtered off nonbiallelic calls 


## Script Usage 


```bash
Usage: ./filter_biallelic_calls.R [options]


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
./filter_biallelic_calls.R --amino_acid_calls ../../data/example2_amino_acid_calls.tsv --out biallelic_example2_amino_acid_calls.tsv.tsv.gz --overwrite

./filter_biallelic_calls.R --amino_acid_calls ../../data/example2_amino_acid_calls.tsv --out biallelic_example2_amino_acid_calls.tsv.tsv.gz --overwrite --out_nonbiallelic nonbiallelic_example2_amino_acid_calls.tsv.tsv.gz
```

