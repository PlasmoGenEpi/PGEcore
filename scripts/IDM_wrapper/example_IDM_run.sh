#! /usr/bin/env bash

# start from the root folder
cd $(git rev-parse --show-toplevel)
cd scripts/IDM_wrapper/

Rscript ./IDM_wrapper.R --aa_calls_input ../../data/example_amino_acid_calls.tsv --model IDM --slaf_output out_slaf.tsv --eps_initial 0.1
