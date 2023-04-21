#!/bin/bash

if [[ $# -lt 2 ]]; then
	echo "Requires at least two arguments: show_unfolded_results.sh <input_file> <output_folder>"
	exit 1
fi

. ./setup.sh
python3 scripts/plotting/unfolding_plots.py $1 -o $2