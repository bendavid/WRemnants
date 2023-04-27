#!/bin/bash

if [[ $# -lt 2 ]]; then
	echo "Requires at least two arguments: print_impacts.sh <input_file> <output_file>"
	exit 1
fi

. ./setup.sh
python3 scripts/combine/printImpacts.py -f $1 
# TODO: Add option to do this in the same call to the function
python3 scripts/combine/pullsAndImpacts.py --oneSidedImpacts -g -f $1 output -o ${2/.html/_grouped.html} --otherExtensions pdf png 
python3 scripts/combine/pullsAndImpacts.py --oneSidedImpacts -f $1 output -o $2 --otherExtensions pdf png -n 50
# Don't limit the number for the html plot
python3 scripts/combine/pullsAndImpacts.py --oneSidedImpacts -f $1 output -o $2
