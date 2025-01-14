#! /usr/bin/env bash

# start from the root folder
cd $(git rev-parse --show-toplevel)
cd scripts/IDM_wrapper/

Rscript ./IDM_wrapper.R --aa_calls ../../data/example_amino_acid_calls.tsv --model IDM --slaf out_slaf.tsv --eps_initial 0.1
