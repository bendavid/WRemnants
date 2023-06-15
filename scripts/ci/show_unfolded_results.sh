#!/bin/bash

if [[ $# -lt 3 ]]; then
	echo "Requires at least three arguments: show_unfolded_results.sh <input_file> <output_path> <output_folder>"
	exit 1
fi

. ./setup.sh
python3 scripts/plotting/unfolding_xsec.py $1 --asimov $1 -o $2 -f unfolding #--eoscp
python3 scripts/plotting/unfolding_plots.py $1 -o $2 -f unfolding #--eoscp
