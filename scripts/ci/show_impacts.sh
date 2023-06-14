#!/bin/bash

if [[ $# -lt 2 ]]; then
	echo "Requires at least two arguments: print_impacts.sh <input_file> <output_file>"
	exit 1
fi

. ./setup.sh
python3 scripts/combine/printImpacts.py -f $1 
<<<<<<< HEAD
# TODO: Add option to do this in the same call to the function
python3 scripts/combine/pullsAndImpacts.py --oneSidedImpacts -g -f $1 output --outFolder $2 -o ${3/.html/_grouped.html} --otherExtensions pdf png --eoscp
python3 scripts/combine/pullsAndImpacts.py --oneSidedImpacts -f $1 output --outFolder $2 -o $3 --otherExtensions pdf png -n 50 --eoscp
=======
python3 scripts/combine/pullsAndImpacts.py --oneSidedImpacts -f $1 output --outFolder $2 -o $3 --otherExtensions pdf png -n 50
>>>>>>> 71c2c47 (Update CI to be consistent with impact plot changes)
# Don't limit the number for the html plot
python3 scripts/combine/pullsAndImpacts.py --oneSidedImpacts -f $1 output --outFolder $2 -o $3 --eoscp
