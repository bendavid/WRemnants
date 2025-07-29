#!/usr/bin/env python3

import argparse

from wremnants.combination_tools import makeTxtresults

parser = argparse.ArgumentParser()


parser.add_argument("-i", "--inputDir", default="/scratch/submit/cms/jbendavi/wmassdev87/data/CMS_Preliminary_2024")
parser.add_argument("-p", "--pdfs", nargs='+', default=["ct18", "ct18z", "msht20", "msht20an3lo", "nnpdf31", "nnpdf40", "pdf4lhc21"])

args = parser.parse_args()

nmap = {}

nmap["ct18"] = "14000_CT18NNLO"
nmap["ct18z"] = "14100_CT18ZNNLO"
nmap["msht20"] = "27400_MSHT20nnlo_as118"
nmap["msht20an3lo"] = "29100_MSHT20an3lo_as118"
nmap["nnpdf31"] = "306000_NNPDF31_nnlo_hessian_pdfas"
nmap["nnpdf40"] = "331600_NNPDF40_nnlo_hessian_pdfas"
nmap["pdf4lhc21"] = "93300_PDF4LHC21_40_pdfas"

smap = {}
smap["ct18"] = "1p0"
smap["ct18z"] = "1p0"
smap["msht20"] = "1p5"
smap["msht20an3lo"] = "1p5"
smap["nnpdf31"] = "3p0"
smap["nnpdf40"] = "5p0"
smap["pdf4lhc21"] = "1p0"

tags = ["pdf", "pdf_infl1"]

for tag in tags:
    for pdf in args.pdfs:
        infile = f"{args.inputDir}/{tag}/{pdf}/WMass_eta_pt_charge/fitresults_123456789_unblind.hdf5"
        asyms = [False, True] if pdf in ["ct18", "ct18z", "msht20", "msht20an3lo"] else [False]
        for asym in asyms:
            tagname = "unscaled" if tag=="pdf_infl1" else f"scale_{smap[pdf]}"
            outfile = f"{nmap[pdf]}_{tagname}"
            if asym:
                outfile += "_asymquad"
            outfile += "_CMS_Preliminary_2024"

            print("in", infile)
            print("out", outfile)

            makeTxtresults(inputFile=infile, outputFile=outfile, asym=asym)

