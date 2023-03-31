#!/bin/bash

if [[ $# -lt 2 ]]; then
	echo "Requires at least two arguments: print_impacts.sh <input_file> <output_file>"
	exit 1
fi

. ./setup.sh
python3 scripts/combine/printImpacts.py -f $1 
# TODO: Add option to do this in the same call to the function
python3 scripts/combine/pullsAndImpacts.py -g -f $1 output -o ${2/.html/_grouped.html}
python3 scripts/combine/pullsAndImpacts.py -f $1 output -o $2
