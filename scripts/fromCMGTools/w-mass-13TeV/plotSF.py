#!/usr/bin/env python3

# plot all SF from the input root file, and make some products for faster usage
# for isolation and antiisolation, new SF can be defined by varying only the data efficiecy or MC one, which is needed to properly account for the anticorrelation between isolated and antiisolated leptons (forcefully anticorrelating the global sf uncertainty is not correct, the 100% anticorrelation holds only for the efficiencies) 
# actually, we verified that the correlation is preserved also for the SF, since isoSF and antiisoSF are found to be generally anticorrelated by more than 95% (even 99% for |eta| < 2.0) for pt < 60 GeV
# no need to do the same also for iso with no trigger, because this splitting is relevant when there are fakes, and for the 2 lepton phase space that background is negligible

# python w-mass-13TeV/plotSF.py testMuonSF/2021-05-31_allSFs_nodz_dxybs.root plots/testNanoAOD/testSF/SFeta0p1_31May2021_nodz_dxybs/globalAndPerEra/ -e 'BtoF,GtoH,B,C,D,E,F,G,H' -n 'trigger,idip,iso,antiiso,isonotrig,antiisonotrig'

import re
import os, os.path
import logging
import argparse
import shutil

## safe batch mode
import sys
args = sys.argv[:]
sys.argv = ['-b']
import ROOT
sys.argv = args
ROOT.gROOT.SetBatch(True)
ROOT.PyConfig.IgnoreCommandLineOptions = True

from copy import *

sys.path.append(os.getcwd() + "/plotUtils/")
from utility import *

logging.basicConfig(level=logging.INFO)

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("rootfile", type=str, nargs=1)
    parser.add_argument("outdir",   type=str, nargs=1, help="Output folder (subfolder 'allSF/' created automatically inside)")
    parser.add_argument("-e", "--era",    type=str, default="BtoF,GtoH,B,C,D,E,F,G,H", help="Comma separated list of eras for SF in histogram name; default: %(default)s")
    parser.add_argument("-n", "--sfnames", type=str, default="trigger,idip,iso,antiiso,isonotrig,antiisonotrig,tracking,reco", help="Comma separated list of SF names inside root file, which will be plotted (trigger uses both plus and minus automatically); default: %(default)s")
    parser.add_argument("--sf-version", dest="sfversions", type=str, default="nominal,dataAltSig", help="SF versions to plot and to use for the products, usually one would use just nominal and dataAltSig to define the systematic variation; default: %(default)s")
    parser.add_argument("--eff-version", dest="effversions", type=str, default="nominal,altSig", help="Efficiency versions to plot (nominal actually has no keyword); default: %(default)s")
    parser.add_argument(     '--nContours', dest='nContours',    default=51, type=int, help='Number of contours in palette. Default is 51')
    parser.add_argument(     '--palette'  , dest='palette',      default=87, type=int, help='Set palette: default is a built-in one, 55 is kRainbow')
    parser.add_argument(     '--invertPalette', dest='invertePalette', action='store_true',   help='Inverte color ordering in palette')
    parser.add_argument(     '--skip-eff', dest='skipEfficiency', action='store_true',   help='Do not plot efficiencies to save time')
    args = parser.parse_args()

    ROOT.TH1.SetDefaultSumw2()

    # products of scale factors
    productsToMake = {"isoTrigPlus"       : ["iso",           "triggerplus",  "idip"], # "tracking"],
                      "isoTrigMinus"      : ["iso",           "triggerminus", "idip"], # "tracking"],
                      "isoNotrig"         : ["isonotrig",                     "idip"], # "tracking"],
                      "noisoTrigPlus"     : [                 "triggerplus",  "idip"], # "tracking"],
                      "noisoTrigMinus"    : [                 "triggerminus", "idip"], # "tracking"],
                      "noisoNotrig"       : [                                 "idip"], # "tracking"],
                      "antiisoTrigPlus"   : ["antiiso",       "triggerplus",  "idip"], # "tracking"],
                      "antiisoTrigMinus"  : ["antiiso",       "triggerminus", "idip"], # "tracking"],
                      "antiisoNotrig"     : ["antiisonotrig",                 "idip"], # "tracking"],
                      "isoOnly"           : ["iso"],
                      "isoNotrigOnly"     : ["isonotrig"],
                      "antiisoOnly"       : ["antiiso"],
                      "antiisoNotrigOnly" : ["antiisonotrig"],
                      "reco"              : ["reco"],
                      "tracking"          : ["tracking"],
                      #"trackingReco"      : ["tracking", "reco"],
                      #"idipTrackingReco"  : ["idip", "tracking", "reco"],
                      #"trigPlusIdipTrackingReco"     : ["triggerplus", "idip", "tracking", "reco"],
                      #"isoTrigPlusIdipTrackingReco"  : ["iso", "triggerplus", "idip", "tracking", "reco"],
                      #"trigMinusIdipTrackingReco"    : ["triggerminus", "idip", "tracking", "reco"],
                      #"isoTrigMinusIdipTrackingReco" : ["iso", "triggerminus", "idip", "tracking", "reco"],
                      #"isoNotrigIdipTrackingReco"    : ["isonotrig", "idip", "tracking", "reco"],
    }

    #productsToMake = {"trackingReco"      : ["reco", "tracking"], # "tracking"],
    #}
    
    eras = args.era.split(',')
    fname = args.rootfile[0]
    outdirOriginal = args.outdir[0] # to keep it below
    if not outdirOriginal.endswith('/'):
        outdirOriginal = outdirOriginal + '/'
    outdirOriginal += "allSF/"
    productSubfolder = "productSF/"
    
    # make folder structure
    foldersToCreate = []
    foldersToCreate.append(outdirOriginal + productSubfolder)
    for era in eras:
        foldersToCreate.append(outdirOriginal + era + "/")
        foldersToCreate.append(outdirOriginal + productSubfolder + era + "/")

    for folder in foldersToCreate:
        createPlotDirAndCopyPhp(folder)
        
    print()
        
    histsSF = {}
    histsEff = {}
    for era in eras:
        histsSF[era] = {}
        histsEff[era] = {"Data": {},
                         "MC"  : {}}

    sf_version = [str(x) for x in args.sfversions.split(",")]
    eff_version = [str(x) for x in args.effversions.split(",")]
    
    names = args.sfnames.split(',')
    f = safeOpenFile(fname)
    for n in names:
        print("")
        print(f"Working point: {n}")
        print("-"*30)
        charges = ["plus", "minus"] if n == "trigger" else ["both"]
        for ch in charges:
            for era in eras:
                tmpEra = "F_preVFP" if era == "F" else era
                for v in sf_version:
                    hname = f"SF2D_{v}_{n}_{tmpEra}_{ch}"
                    hkey = f"{v}_{n}" 
                    hkey += (ch if n == "trigger" else "")
                    print(f"   {hkey} -> {hname}")
                    histsSF[era][hkey] = safeGetObject(f, hname)
                for v in eff_version:
                    for dataMC in ["Data", "MC"]:
                        realv = "" if v == "nominal" else f"{v}_"
                        hname = f"eff{dataMC}_{realv}{n}_{tmpEra}_{ch}"
                        hkey = f"{v}_{n}" 
                        hkey += (ch if n == "trigger" else "")
                        print(f"   {hkey} -> {hname}")
                        histsEff[era][dataMC][hkey] = safeGetObject(f, hname)
    f.Close()

    canvas = ROOT.TCanvas("canvas","",800,800)

    canvas_unroll = ROOT.TCanvas("canvas_unroll","",3000,800) 
    leftMargin = 0.06
    rightMargin = 0.01
    bottomMargin = 0.12
    canvas_unroll.SetTickx(1)
    canvas_unroll.SetTicky(1)
    canvas_unroll.cd()
    canvas_unroll.SetLeftMargin(leftMargin)
    canvas_unroll.SetRightMargin(rightMargin)
    canvas_unroll.cd()
    canvas_unroll.SetBottomMargin(bottomMargin)

    # first plot efficiencies
    if not args.skipEfficiency:
        for era in eras:
            outdir = outdirOriginal + era + "/"
            for dataMC in ["Data", "MC"]:
                for n in list(histsEff[era][dataMC].keys()):
                    drawCorrelationPlot(histsEff[era][dataMC][n], "muon #eta", "muon p_{T} (GeV)", f"{dataMC} efficiency",
                                        f"muonEff{dataMC}_{n}", plotLabel="ForceTitle", outdir=outdir,
                                        smoothPlot=False, drawProfileX=False, scaleToUnitArea=False,
                                        draw_both0_noLog1_onlyLog2=1, passCanvas=canvas,
                                        nContours=args.nContours, palette=args.palette, invertePalette=args.invertePalette)
                    # abs. uncertainty (only on nominal, it should be the same for all histograms, hoping the SF and efficiencies were sane in this configuration)
                    if "nominal" in n:
                        drawCorrelationPlot(histsEff[era][dataMC][n], "muon #eta", "muon p_{T} (GeV)", f"Abs. uncertainty on {dataMC} eff",
                                            f"absUnc_muonEff{dataMC}_{n}", plotLabel="ForceTitle", outdir=outdir+"absoluteStatUncertainty/",
                                            smoothPlot=False, drawProfileX=False, scaleToUnitArea=False,
                                            draw_both0_noLog1_onlyLog2=1, passCanvas=canvas, plotError=True,
                                            nContours=args.nContours, palette=args.palette, invertePalette=args.invertePalette)
                        ## rel. uncertainty
                        drawCorrelationPlot(histsEff[era][dataMC][n], "muon #eta", "muon p_{T} (GeV)", f"Rel. uncertainty on {dataMC} eff",
                                            f"relUnc_muonEff{dataMC}_{n}", plotLabel="ForceTitle", outdir=outdir+"relativeStatUncertainty/",
                                            smoothPlot=False, drawProfileX=False, scaleToUnitArea=False,
                                            draw_both0_noLog1_onlyLog2=1, passCanvas=canvas, plotRelativeError=True,
                                            nContours=args.nContours, palette=args.palette, invertePalette=args.invertePalette)

                        
    prodHistsSF = {}
    for era in eras:
        prodHistsSF[era] = {}

    namesForCheck = [x for x in names if x != "trigger"]
    if "trigger" in names:
        namesForCheck += ["triggerplus", "triggerminus"]
        
    for key in list(productsToMake.keys()):
        if not all(x in namesForCheck for x in productsToMake[key]):
            print()
            missingFactors = ",".join([x for x in productsToMake[key] if x not in namesForCheck])
            print(f"Warning: skipping product {key} because these factors are missing: {missingFactors}")
            print()
            continue
        for era in eras:
            for sfv in sf_version:
                prodname = f"fullSF2D_{sfv}_{key}_{era}"
                for i,basename in enumerate(productsToMake[key]): 
                    name = f"{sfv}_{basename}"
                    if i == 0:
                        stringProduct = basename
                        prodHistsSF[era][prodname] = copy.deepcopy(histsSF[era][name].Clone(prodname))
                    else:
                        stringProduct = stringProduct + "*" + basename
                        if histsSF[era][name].GetNbinsY() == 1:
                            multiplyByHistoWith1ptBin(prodHistsSF[era][prodname], histsSF[era][name])
                        elif not prodHistsSF[era][prodname].Multiply(histsSF[era][name]):
                            print(f"ERROR in multiplication for prodHistsSF[{era}][{prodname}] with {name}")
                            print(f"Nbins(X, Y) = {histsSF[era][name].GetNbinsX()},{histsSF[era][name].GetNbinsY()} ")
                            quit()
                prodHistsSF[era][prodname].SetTitle(f"{stringProduct}")            
                print(f"{era}: {sfv} -> {stringProduct}")
                
    fullout = outdirOriginal + "scaleFactorProduct.root" 
    f = ROOT.TFile.Open(fullout, "RECREATE")
    # put new ones in this files, without copying original histograms (although some are actually bare copies when the product has a single factor)
    if not f or not f.IsOpen():
        raise RuntimeError(f"Error when opening file {fullout}")
    f.cd()
    for era in eras:
        # save the product of scale factors
        for n in list(prodHistsSF[era].keys()):
            prodHistsSF[era][n].Write(n)
    f.Close()

    minmax = {"triggerplus"  : "0.65,1.15",
              "triggerminus" : "0.65,1.15",
              "idip"         : "0.95,1.01",
              "iso"          : "0.975,1.025",
              "antiiso"      : "0.7,1.25",
              "isonotrig"    : "0.97,1.03",
              "antiisonotrig": "0.7,1.25",
              "tracking"     : "0.98,1.01",
              "reco"         : "0.94,1.02",
    }
    
    for era in eras:

        outdir = outdirOriginal + era + "/"

        for n in list(histsSF[era].keys()):
            version = str(n.split("_")[0])
            ntmp    = str(n.split("_")[1])
            if ntmp in minmax.keys():
                zrange = f"::{minmax[ntmp]}"
            else:
                zrange = ""
                        
            drawCorrelationPlot(histsSF[era][n], "muon #eta", "muon p_{T} (GeV)", f"Data/MC scale factor{zrange}",
                                f"muonSF_{n}", plotLabel="ForceTitle", outdir=outdir,
                                smoothPlot=False, drawProfileX=False, scaleToUnitArea=False,
                                draw_both0_noLog1_onlyLog2=1, passCanvas=canvas,
                                nContours=args.nContours, palette=args.palette, invertePalette=args.invertePalette)
            # abs. uncertainty (only on nominal, it should be the same for all histograms, hoping the SF and efficiencies were sane in this configuration)
            if "nominal" in n:
                drawCorrelationPlot(histsSF[era][n], "muon #eta", "muon p_{T} (GeV)", f"Abs. uncertainty on SF",
                                    f"absUnc_muonSF_{n}", plotLabel="ForceTitle", outdir=outdir+"absoluteStatUncertainty/",
                                    smoothPlot=False, drawProfileX=False, scaleToUnitArea=False,
                                    draw_both0_noLog1_onlyLog2=1, passCanvas=canvas, plotError=True,
                                    nContours=args.nContours, palette=args.palette, invertePalette=args.invertePalette)
                ## rel. uncertainty
                drawCorrelationPlot(histsSF[era][n], "muon #eta", "muon p_{T} (GeV)", f"Rel. uncertainty on SF",
                                    f"relUnc_muonSF_{n}", plotLabel="ForceTitle", outdir=outdir+"relativeStatUncertainty/",
                                    smoothPlot=False, drawProfileX=False, scaleToUnitArea=False,
                                    draw_both0_noLog1_onlyLog2=1, passCanvas=canvas, plotRelativeError=True,
                                    nContours=args.nContours, palette=args.palette, invertePalette=args.invertePalette)
                # plot systs
                hnomi = copy.deepcopy(histsSF[era][n].Clone("hnomi"))
                unrolledSF = [unroll2Dto1D(histsSF[era][n], newname=f"unrolled_{n}", cropNegativeBins=False)]
                unrolledLeg = [f"Nominal {ntmp} scale factors"]
                for n2 in list(histsSF[era].keys()):
                    if "nominal" in n2:
                        continue
                    v2 = n2.split("_")[0]
                    n2tmp = n2.split("_")[1]
                    if ntmp != n2tmp:
                        continue
                    hsyst = copy.deepcopy(histsSF[era][n2].Clone(f"hsyst_{n2}"))
                    unrolledSF.append(unroll2Dto1D(hsyst, newname=f"unrolled_{n2}", cropNegativeBins=False))
                    unrolledLeg.append(f"{v2} syst")
                    hsyst.Add(hnomi, -1)
                    drawCorrelationPlot(hsyst, "muon #eta", "muon p_{T} (GeV)", f"{v2} syst unc. on SF",
                                        f"muonSF_{n2}_syst", plotLabel="ForceTitle", outdir=outdir+"absoluteSystUncertainty/",
                                        smoothPlot=False, drawProfileX=False, scaleToUnitArea=False,
                                        draw_both0_noLog1_onlyLog2=1, passCanvas=canvas,
                                        nContours=args.nContours, palette=args.palette, invertePalette=args.invertePalette)
                ptBinRanges = []
                for ipt in range(hnomi.GetNbinsY()):
                    ptBinRanges.append("#splitline{{[{ptmin},{ptmax}]}}{{GeV}}".format(ptmin=int(hnomi.GetYaxis().GetBinLowEdge(ipt+1)),
                                                                                       ptmax=int(hnomi.GetYaxis().GetBinLowEdge(ipt+2))))

                drawNTH1(unrolledSF, unrolledLeg, "Unrolled tag-and-probe eta-p_{T} bin", f"Scale factor", f"scaleFactorAndUncertainty_{ntmp}", outdir+"unrolled/",
                         leftMargin=0.06, rightMargin=0.01, labelRatioTmp="Syst/Nominal",
                         legendCoords="0.2,0.8,0.9,0.98;5", lowerPanelHeight=0.5, skipLumi=True, passCanvas=canvas_unroll,
                         drawVertLines="{a},{b}".format(a=hnomi.GetNbinsY(),b=hnomi.GetNbinsX()), yAxisExtendConstant=1.4,
                         textForLines=ptBinRanges, transparentLegend=False, drawErrorAll=False,
                         onlyLineColor=True, useLineFirstHistogram=True, setRatioRangeFromHisto=True, setOnlyLineRatio=False,
                         lineWidth=1, ytextOffsetFromTop=0.2, useMultiHistRatioOption=True)

                
        outdir = outdirOriginal + productSubfolder + era + "/"

        for n in list(prodHistsSF[era].keys()):
            drawCorrelationPlot(prodHistsSF[era][n], "muon #eta", "muon p_{T} (GeV)", "Data/MC scale factor product",
                                f"{n}", plotLabel="ForceTitle", outdir=outdir,
                                smoothPlot=False, drawProfileX=False, scaleToUnitArea=False,
                                draw_both0_noLog1_onlyLog2=1, passCanvas=canvas,
                                nContours=args.nContours, palette=args.palette, invertePalette=args.invertePalette)
            # plot absolute error
            drawCorrelationPlot(prodHistsSF[era][n], "muon #eta", "muon p_{T} (GeV)", "Abs. uncertainty on SF product",
                                f"absUnc_{n}", plotLabel="ForceTitle", outdir=outdir+"absoluteStatUncertainty/",
                                smoothPlot=False, drawProfileX=False, scaleToUnitArea=False,
                                draw_both0_noLog1_onlyLog2=1, passCanvas=canvas, plotError=True,
                                nContours=args.nContours, palette=args.palette, invertePalette=args.invertePalette)
            ## plot relative error
            drawCorrelationPlot(prodHistsSF[era][n], "muon #eta", "muon p_{T} (GeV)", "Rel. uncertainty on SF product",
                                f"relUnc_{n}", plotLabel="ForceTitle", outdir=outdir+"relativeStatUncertainty/",
                                smoothPlot=False, drawProfileX=False, scaleToUnitArea=False,
                                draw_both0_noLog1_onlyLog2=1, passCanvas=canvas, plotRelativeError=True,
                                nContours=args.nContours, palette=args.palette, invertePalette=args.invertePalette)
