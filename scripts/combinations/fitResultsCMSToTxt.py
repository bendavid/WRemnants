#!/usr/bin/env python3

import argparse

from wremnants.combination_tools import makeTxtresults

parser = argparse.ArgumentParser()

parser.add_argument("-i", "--inputFile")
parser.add_argument("-o", "--outputFile")

args = parser.parse_args()

makeTxtresults(inputFile = args.inputFile, outputFile = args.outputFile)
