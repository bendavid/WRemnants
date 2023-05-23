## example
# python tests/testShapesIsoMtRegions.py path/to/pickle/file.pkl.lz4 /path/to/plot/folder/testShapesIsoMtRegions/ [-c minus] [--isoMtRegion 2]

from wremnants.datasets.datagroups2016 import make_datagroups_2016
from wremnants import histselections as sel
#from wremnants import plot_tools,theory_tools,syst_tools
from utilities import boostHistHelpers as hh
from utilities import common, logging

import narf
import wremnants
from wremnants import theory_tools,syst_tools,theory_corrections
import hist

import numpy as np
from utilities import input_tools

import lz4.frame

import argparse
import os
import shutil
import re

## safe batch mode
import sys
args = sys.argv[:]
sys.argv = ['-b']
import ROOT
sys.argv = args
ROOT.gROOT.SetBatch(True)
ROOT.PyConfig.IgnoreCommandLineOptions = True

from copy import *

#sys.path.append(os.getcwd() + "/plotUtils/")
#from utility import *
from scripts.analysisTools.plotUtils.utility import *
from scripts.analysisTools.w_mass_13TeV.plotPrefitTemplatesWRemnants import plotPrefitHistograms

sys.path.append(os.getcwd())

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("inputfile", type=str, nargs=1, help="Input file with histograms (pkl.lz4 or hdf5 file)")
    parser.add_argument("outdir",   type=str, nargs=1, help="Output folder")
    parser.add_argument("-v", "--verbose", type=int, default=3, choices=[0,1,2,3,4], help="Set verbosity level with logging, the larger the more verbose");
    parser.add_argument("-n", "--baseName", type=str, help="Histogram name in the file (e.g., 'nominal')", default="nominal")
    parser.add_argument(     '--nContours', dest='nContours',    default=51, type=int, help='Number of contours in palette. Default is 51')
    parser.add_argument(     '--palette'  , dest='palette',      default=87, type=int, help='Set palette: default is a built-in one, 55 is kRainbow')
    parser.add_argument(     '--invertPalette', dest='invertePalette', action='store_true',   help='Inverte color ordering in palette')
    parser.add_argument('-p','--processes', default=None, nargs='*', type=str,
                        help='Choose what processes to plot, otherwise all are done')
    parser.add_argument('-c','--charges', default=["plus", "minus"], choices=["plus", "minus"], nargs='+', type=str,
                        help='Choose what charge to plot')
    parser.add_argument("--isoMtRegion", type=int, nargs='+', default=[0,1,2,3], choices=[0,1,2,3], help="Integer index for iso-Mt regions to plot (conversion is index = passIso * 1 + passMT * 2 as in common.getIsoMtRegionFromID)");
    args = parser.parse_args()

    logger = logging.setup_logger(os.path.basename(__file__), args.verbose, True)
    # if 0:
    #     logger.critical("TEST LOGGER CRITICAL")
    #     logger.error("TEST LOGGER ERROR")
    #     logger.warning("TEST LOGGER WARNING")
    #     logger.info("TEST LOGGER INFO")
    #     logger.debug("TEST LOGGER DEBUG")
    #     quit()

    fname = args.inputfile[0]
    outdir = args.outdir[0]
    createPlotDirAndCopyPhp(outdir)
        
    ROOT.TH1.SetDefaultSumw2()

    canvas = ROOT.TCanvas("canvas", "", 800, 700)
    cwide = ROOT.TCanvas("cwide","",2400,600)                      
    adjustSettings_CMS_lumi()
    canvas1D = ROOT.TCanvas("canvas1D", "", 800, 900)

    xAxisName = "Muon #eta"
    yAxisName = "Muon p_{T} (GeV)"

    for isoMtID in args.isoMtRegion:

        ret = common.getIsoMtRegionFromID(isoMtID)
        passMT = ret["passMT"]
        passIso = ret["passIso"]
        passMT_str = "passMT" if passMT else "failMT"
        passIso_str = "passIso" if passIso else "failIso"
        selectOp = None # currently same for all processes
        if   not passIso and not passMT:
            selectOp = sel.histWmass_failMT_failIso
        elif     passIso and not passMT:
            selectOp = sel.histWmass_failMT_passIso
        elif not passIso and     passMT:
            selectOp = sel.histWmass_passMT_failIso
        else:
            selectOp = sel.histWmass_passMT_passIso
            # customized for fakes later on
            
        groups = make_datagroups_2016(fname)
        datasets = groups.getNames()
        if args.processes is not None and len(args.processes):
            datasets = list(filter(lambda x: x in args.processes, datasets))
        logger.info(f"Will plot datasets {datasets}")

        nominalName = args.baseName.rsplit("_", 1)[0]
        groups.setNominalName(nominalName)
        # set region
        groups.setSelectOp(selectOp)
        if passIso and passMT and "Fake" in datasets:
            groups.setSelectOp(sel.fakeHistABCD, processes=["Fake"])
        groups.loadHistsForDatagroups(nominalName, syst="", procsToRead=datasets)

        histInfo = groups.getDatagroups()
        #print(histInfo)
        print("-"*30)
        print(f"Region ID = {isoMtID}")
        print("-"*30)
        rootHists = {}

        charges = args.charges

        hist2D = {c : {} for c in charges}
        for d in datasets:
            hnarf = histInfo[d].hists[args.baseName]
            #print(f"{d}: {hnarf.sum()}")
            rootHists[d] = narf.hist_to_root(hnarf) # this is a TH3 with eta, pt, charge
            for charge in charges:
                rootHists[d].SetName(f"nominal_isoMT_{isoMtID}_{d}")
                chargeBin = 1 if charge == "minus" else 2
                h = getTH2fromTH3(rootHists[d], f"{rootHists[d].GetName()}_{charge}", chargeBin, chargeBin)
                h.SetTitle(f"{d} {charge}: {passMT_str} {passIso_str}")
                drawCorrelationPlot(h, xAxisName, xAxisName, f"Events",
                                    f"{h.GetName()}", plotLabel="ForceTitle", outdir=outdir,
                                    smoothPlot=False, drawProfileX=False, scaleToUnitArea=False,
                                    draw_both0_noLog1_onlyLog2=1, passCanvas=canvas,
                                    nContours=args.nContours, palette=args.palette, invertePalette=args.invertePalette)
                hist2D[charge][f"{d}_{charge}"] = h
                print()

        if "Data" in datasets:
            outdir_isoMtID = f"{outdir}/isoMtID_{isoMtID}/"
            createPlotDirAndCopyPhp(outdir_isoMtID)
            for c in charges:
                hdata2D = hist2D[c][f"Data_{c}"]        
                hmc2D = [hist2D[c][key] for key in hist2D[c].keys() if key.split("_")[0] not in ["Data", "QCD"]]
                outdir_dataMC = f"{outdir_isoMtID}dataMC_{c}/"
                createPlotDirAndCopyPhp(outdir_dataMC)
                plotPrefitHistograms(hist2D[c][f"Data_{c}"], hmc2D, outdir_dataMC, xAxisName=xAxisName, yAxisName=yAxisName,
                                     chargeLabel=c, canvas=canvas, canvasWide=cwide, canvas1D=canvas1D)
                
    # processes = ["WplusmunuPostVFP"]
    # hnomi = {}
    # for p in processes:
    #     hnomi[p] = input_tools.read_and_scale(fname, p, args.baseName, True)
    #     print(hnomi[p].axes)
    # quit()
