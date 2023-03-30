#!/bin/bash

if [[ $# -lt 2 ]]; then
	echo "Requires at least two arguments: print_impacts.sh <input_file> <output_file>"
	exit 1
fi

. ./setup.sh
python3 scripts/combine/printImpacts.py -f $1 
python3 scripts/combine/pullsAndImpacts.py -g -f $1 output -o $2
