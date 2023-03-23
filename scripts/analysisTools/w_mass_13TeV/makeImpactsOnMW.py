#!/bin/env python

# plot impacts on mW or mZ from groups of nuisance parameters

# example for wlike (need --wlike)
# python w-mass-13TeV/makeImpactsOnMW.py /scratch/mciprian/CombineStudies/Wlike/20Sept2022/qcdScale_byHelicityPt/nominal/fit/hessian/fitresults_123456789_Asimov_bbb1_cxs0.root  -o plots/fromMyWremnants/Wlike_fit/10Sept2022/qcdScale_byHelicityPt/makeImpactsOnMW/ --set-stat 0.071 --showTotal --scaleToMeV --wlike
#
# wlike with charge + (assuming root file with fit results exists)
# python w-mass-13TeV/makeImpactsOnMW.py /scratch/mciprian/CombineStudies/Wlike/20Sept2022/qcdScale_byHelicityPt/nominal/fit/hessian/fitresults_123456789_Asimov_onlyplus_bbb1_cxs0.root  -o plots/fromMyWremnants/Wlike_fit/10Sept2022/qcdScale_byHelicityPt/makeImpactsOnMW/ --set-stat 0.100 --showTotal --scaleToMeV --postfix chargePlus --wlike


import os, datetime, re, operator, math
import argparse
import shutil

from array import array
from copy import *

## safe batch mode
import sys
args = sys.argv[:]
sys.argv = ['-b']
import ROOT
sys.argv = args
ROOT.gROOT.SetBatch(True)
ROOT.PyConfig.IgnoreCommandLineOptions = True

import utilitiesCMG
utilities = utilitiesCMG.util()

from scripts.analysisTools.plotUtils.utility import *

def readNuisances(args, infile=None):

    if infile is None:
        infile = args.rootfile[0]

    print(f"Starting with file {infile} ...")

    nuis = []
    if args.nuisgroups == "ALL": 
        nuis = ["ALL"]
    else:
        nuis = list(args.nuisgroups.split(','))

    massNuisanceName = "massShift{s}MeV".format(s=int(args.prefitUncertainty))
    valuesAndErrors = utilities.getFromHessian(infile,params=[massNuisanceName])
    totalUncertainty = valuesAndErrors[massNuisanceName][1] - valuesAndErrors[massNuisanceName][0]
    if args.scaleToMeV:
        totalUncertainty *= args.prefitUncertainty
        print("Total m%s uncertainty: %2.1f MeV" % (boson, totalUncertainty))
    else:
        print("Total m%s uncertainty: %2.3f (\% of prefit)" % (boson, totalUncertainty))
        
    group = 'group_' if len(args.nuisgroups) else ''
    th2name = 'nuisance_{group}impact_nois'.format(group=group)
    hessfile = ROOT.TFile(infile,'read')
    impMat = hessfile.Get(th2name)
    if impMat==None:
        print("ERROR: Cannot find the impact TH2 named ",th2name," in the input file. Maybe you didn't run --doImpacts?\nSkipping.")
        sys.exit()

    if args.excludeNuisgroups:
        matchExclude = re.compile(args.excludeNuisgroups)
        
    print("Histograms loaded successfully ...")
    nuisGroup_nameVal = {}
    for iy in range(1,impMat.GetNbinsY()+1):
        label = impMat.GetYaxis().GetBinLabel(iy)
        if args.excludeNuisgroups and matchExclude.match(label):
            continue
        if nuis[0] != "ALL":
            if label in nuis:
                nuisGroup_nameVal[label] = impMat.GetBinContent(1,iy)
        else:
            nuisGroup_nameVal[label] = impMat.GetBinContent(1,iy)
    return totalUncertainty,nuisGroup_nameVal


if __name__ == "__main__":

    ROOT.gStyle.SetOptStat(0)

    parser = argparse.ArgumentParser()
    parser.add_argument("rootfile", type=str, nargs=1)
    parser.add_argument('-o','--outdir',     dest='outdir',     default='',   type=str, help='output directory to save the matrix')    
    parser.add_argument(     '--nuisgroups', dest='nuisgroups', default='ALL',   type=str, help='nuis groups for which you want to show the impacts (can pass comma-separated list to make all of them one after the other). Use full name, no regular expressions. By default, all are made')
    parser.add_argument('-x', '--excludeNuisgroups', dest='excludeNuisgroups', default=None,   type=str, help='Regular expression for nuisances to be excluded')
    parser.add_argument(     '--set-stat'  , dest='setStat',    default=-1.0, type=float, help='If positive, use this value for stat (this is before scaling to MeV) until combinetf is fixed')
    parser.add_argument(     '--postfix',     dest='postfix',     default='',   type=str, help='postfix for the output name')
    parser.add_argument(     '--canvasSize', dest='canvasSize', default='800,1200', type=str, help='Pass canvas dimensions as "width,height" ')
    # parser.add_argument(     '--draw-option', dest='drawOption', default='COLZ TEXT', type=str, help='Options for drawing TH2')
    parser.add_argument(     '--margin',     dest='margin',     default='', type=str, help='Pass canvas margin as "left,right,top,bottom" ')
    parser.add_argument(     '--nContours', dest='nContours',    default=51, type=int, help='Number of contours in palette. Default is 51 (let it be odd, so the central strip is white if not using --abs-value and the range is symmetric)')
    parser.add_argument(     '--palette'  , dest='palette',      default=0, type=int, help='Set palette: default is a built-in one, 55 is kRainbow')
    parser.add_argument(     '--invertPalette', dest='invertePalette' , default=False , action='store_true',   help='Inverte color ordering in palette')
    parser.add_argument(     '--scaleToMeV', dest='scaleToMeV' , default=False , action='store_true',   help='Report numbers in terms of uncertainty on mW in MeV (default is to report percentage of prefit uncertainty)')
    parser.add_argument(     '--showTotal', dest='showTotal' , default=False , action='store_true',   help='Show total uncertainty in plot')
    parser.add_argument(     '--prefitUncertainty'  , dest='prefitUncertainty',      default=100.0, type=float, help='prefit uncertainty on mW in MeV')
    parser.add_argument(     '--wlike', dest='isWlike', action="store_true", default=False, help="impacts for W-like analysis (it prints mZ accordingly). Default is Wmass");
    parser.add_argument(     '--compareFile', dest='compareFile', default='', type=str, help='Additional file to compare impacts with (must have the same impact labels)')
    parser.add_argument(     '--set-stat-alt'  , dest='setStatAlt',    default=-1.0, type=float, help='If positive, use this value for stat of the file passed with compareFile')
    parser.add_argument(     '--legendEntries', nargs=2, type=str, help="Legend entries when comparing files", default=["Nominal", "Alternate"])
    args = parser.parse_args()

    # palettes:
    # 69 + inverted, using TColor::InvertPalette(), kBeach
    # 70 + inverted, using TColor::InvertPalette(), kBlackBody
    # 109, kCool
    # 57 kBird
    # 55 kRainBow
    # 62 kColorPrintableOnGrey
    # 73 kCMYK
    # 58 kCubehelix
    # 68 kAvocado
    # 111 kGistEarth + inverted
    # 87 kLightTemperature
    # 85 kIsland
    # 56 kInvertedDarkBodyRadiator
    # 100 kSolar + inverted


    # if len(args.nuis) and len(args.nuisgroups):
    #     print('You can specify either single nuis (--nuis) or poi groups (--nuisgroups), not both!')
    #     sys.exit()
    ROOT.TColor.CreateGradientColorTable(3,
                                         array ("d", [0.00, 0.50, 1.00]),
                                         ##array ("d", [1.00, 1.00, 0.00]),
                                         ##array ("d", [0.70, 1.00, 0.34]),
                                         ##array ("d", [0.00, 1.00, 0.82]),
                                         array ("d", [0.00, 1.00, 1.00]),
                                         array ("d", [0.34, 1.00, 0.65]),
                                         array ("d", [0.82, 1.00, 0.00]),
                                         255,  0.95)


    absValue = False
    if len(args.nuisgroups): absValue = True

    if absValue:    
        ROOT.TColor.CreateGradientColorTable(2,
                                             array ("d", [0.00, 1.00]),
                                             ##array ("d", [1.00, 1.00, 0.00]),
                                             ##array ("d", [0.70, 1.00, 0.34]),
                                             ##array ("d", [0.00, 1.00, 0.82]),
                                             array ("d", [1.00, 1.00]),
                                             array ("d", [1.00, 0.65]),
                                             array ("d", [1.00, 0.00]),
                                             255,  0.95)

    boson = "Z" if args.isWlike else "W"

    createPlotDirAndCopyPhp(args.outdir)

    compare = True if len(args.compareFile) else False
        
    totalUncertainty_mW, nuisGroup_nameVal = readNuisances(args, args.rootfile[0])
    if args.setStat > 0.0:
        nuisGroup_nameVal["stat"] = args.setStat

    if compare:
        totalUncertainty_mW_alt, nuisGroup_nameVal_alt = readNuisances(args, args.compareFile)
        if args.setStatAlt > 0.0:
            nuisGroup_nameVal_alt["stat"] = args.setStatAlt

    sortedGroups = sorted(nuisGroup_nameVal.keys(), key= lambda x: nuisGroup_nameVal[x])

    ROOT.gStyle.SetPaintTextFormat('2.1f' if args.scaleToMeV else '0.3f')

    print("Creating output histogram ...")
    # add 1 more bin for total
    nbins = len(sortedGroups)
    if args.showTotal:
        nbins += 1

    # new version
    h1 = ROOT.TH1D("impactsOnMw_chart","impacts of nuisance groups on m_{%s}" % boson,nbins,0,nbins)
    if compare:
        h2 = h1.Clone("impactsOnMw_chart_alt")
    h1.GetYaxis().SetTitle("impacts on m_{{{boson}}} {units}".format(boson=boson, units="[MeV]" if args.scaleToMeV else ""))
    h1.GetYaxis().SetTitleOffset(1.05)
    h1.GetYaxis().SetTitleSize(0.045)
    h1.GetYaxis().SetLabelSize(0.04)
    #h1.GetYaxis().SetTitle("")
    for ik,k in enumerate(sortedGroups):
        bincontent = nuisGroup_nameVal[k] if not args.scaleToMeV else nuisGroup_nameVal[k] * args.prefitUncertainty
        print("%s: %2.1f" % (k,bincontent))
        #print("%s: %2.1f" % (k,bincontent))
        h1.GetXaxis().SetBinLabel(ik+1,k)
        h1.SetBinContent(ik+1,bincontent)
        if compare:
            bincontentAlt = nuisGroup_nameVal_alt[k] if not args.scaleToMeV else nuisGroup_nameVal_alt[k] * args.prefitUncertainty
            h2.SetBinContent(ik+1,bincontentAlt)
    if args.showTotal:
        h1.GetXaxis().SetBinLabel(nbins,"total")
        h1.SetBinContent(nbins,totalUncertainty_mW)
        if compare:
            h2.SetBinContent(nbins,totalUncertainty_mW_alt)
            
    h1.GetXaxis().SetTitleOffset(1.2)
    h1.GetXaxis().SetTitleSize(0.05)
    h1.GetXaxis().SetLabelSize(0.05)

    h1.SetBarWidth(0.8);
    h1.SetBarOffset(0.1);
    h1.SetFillColor(ROOT.kGreen-5)
    h1.SetLineColor(ROOT.kBlack)
    #h1.SetFillStyle(3001)
    h1.SetFillColorAlpha(ROOT.kGreen-5, 0.5)
    if compare:
        h2.SetBarWidth(0.4);
        h2.SetBarOffset(0.2);
        h2.SetFillColor(ROOT.kPink-6)
        h2.SetLineColor(ROOT.kBlack)
        #h2.SetFillStyle(3001)
        h2.SetFillColorAlpha(ROOT.kPink-6, 0.5)
        h1.GetYaxis().SetRangeUser(0.0, 1.1* max(totalUncertainty_mW_alt, totalUncertainty_mW))
        
    cw,ch = args.canvasSize.split(',')
    c1 = ROOT.TCanvas("c1", "", int(cw), int(ch))
    #c1.SetFillColor(42)
    c1.SetGridx()
    c1.SetGridy()
    #ROOT.gPad.SetFrameFillColor(33)
    
    clm = 0.4
    crm = 0.12
    cbm = 0.1
    ctm = 0.1
    if args.margin:
        clm,crm,ctm,cbm = (float(x) for x in args.margin.split(','))
    c1.SetLeftMargin(clm)
    c1.SetRightMargin(crm)
    c1.SetBottomMargin(cbm)
    c1.SetTopMargin(ctm)
    c1.SetTickx(1)
    c1.SetTicky(1)
    h1.Draw("hbar1")
    if compare:
        h2.Draw("hbar1 SAME")
        leg = ROOT.TLegend(0.1,0.9,0.92,0.95)
        leg.SetNColumns(2)
        leg.SetFillColor(0)
        leg.SetFillColorAlpha(0,0.6)
        leg.SetShadowColor(0)
        leg.SetLineColor(0)
        leg.SetBorderSize(0)
        leg.SetTextFont(42)
        leg.SetTextSize(0.035)
        leg.AddEntry(h1, args.legendEntries[0], "LF")
        leg.AddEntry(h2, args.legendEntries[1], "LF")
        leg.Draw("SAME")
        
    #hval = copy.deepcopy(h1.Clone("hval"))
    hval = h1.Clone("hval")
    hval.Reset("ICESM")
    hval.GetXaxis().SetTitleOffset(1.2)
    hval.GetXaxis().SetTitleSize(0.05)
    hval.GetXaxis().SetLabelSize(0.05)
    lat = ROOT.TLatex()
    #lat.SetNDC();
    lat.SetTextFont(42)        
    lat.SetTextSize(0.035)
    c1.Update()
    xtex = 1.05 * (ROOT.gPad.GetUxmax() - ROOT.gPad.GetUxmin())
    step = (ROOT.gPad.GetUymax() - ROOT.gPad.GetUymin()) / h1.GetNbinsX()
    ytex = ROOT.gPad.GetUymin() + step/2
    for i in range(1, 1 + h1.GetNbinsX()):
        hval.GetXaxis().SetBinLabel(i, str(round(h1.GetBinContent(i), 1 if args.scaleToMeV else 3)))
        hval.SetBinContent(i, 0.0)
        lat.DrawLatex(xtex, ytex + step * (i-1), hval.GetXaxis().GetBinLabel(i))
    hval.Draw("AXIS X+ SAME")
        
    postfix = args.postfix
    if len(postfix) and  not postfix.startswith("_"):
        postfix = "_" + postfix
    smallBoson = "z" if args.isWlike else "w"
    for ext in ["pdf", "png"]:
        c1.SaveAs(args.outdir + f"/impactsOnM{smallBoson}{postfix}_chart.{ext}")
