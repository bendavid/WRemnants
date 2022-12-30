#!/bin/env pyth

### Wrapper to run all plotting scripts in one go out of the fit results
### This runs the minimal command, which is usually sufficient
### Each script has additional options for more customization, please check them

import ROOT, os, sys, re, array
sys.path.append(os.getcwd() + "/plotUtils/")
from plotUtils.utility import safeSystem

# general settings
dryRun = 0
isWlike = 1
skipData = 1 # or set fits = ["Asimov"]
onlyData = 0 # or set fits = ["Data"]
fits = ["Asimov", "Data"]

# what to plot
skipHistograms = 0 # prefit histograms, can't be made also before running the fit
skipImpacts = 0
skipNuisances = 0
skipSystRatios = 1
skipPostfitHistograms = 1 # prefit and postfit histograms, from fitresults.root


## These should not be touched
what = "ZMassWLike" if isWlike else "WMass"

# input and output folders
# TODO: need a general way so that every user can use the same logic for paths
mainInputPath = f"/scratch/mciprian/CombineStudies/{what}/smoothSF/muonCorr_none/scetlibCorr_nnpdf31/byHelicityPtCharge/" # contains the TH2 histograms for combine
subFolder = "nominal" # contains the final cards and fit results
mainOutdir = f"/eos/user/m/mciprian/www/WMassAnalysis/fromMyWremnants/{what}_fit/smoothSF/muonCorr_none/scetlibCorr_nnpdf31/byHelicityPtCharge/" # where to store plots
combineInputFile = f"{mainInputPath}/{what}CombineInput.root" 
useSmoothSF = True if "smoothSF" in mainInputPath else False # FIXME: a bit hardcoded

##############################
# utility functions used here
def printText(text):
    print("")
    print("="*30)
    print(text)
    print("")

##############################
# dictionaries to customize some specific plots
#
#
# postfix for plot name and regular expression to pick histogram variations to plot
# there are infinite combinations, so this is mainly an example for Up variations
systRatiosDict = {"effStat_trigger_eta20_plus"  : ".*effStat.*_trigger_eta20.*q1.*Up",
                  "pdfsAndAlphaS"  : ".*pdf(12|25|50|.*AlphaS).*Up",
}

#
# unique string for output plot name and regular expression for nuisances used to plot pulls and constraints
diffNuisanceDict = {"effStat_triggerPlus"  : ".*effStat.*_trigger.*q1",
                    "effStat_triggerMinus" : ".*effStat.*_trigger.*q0",
                    "effStat_idipPlus"  : ".*effStat.*_idip.*q1",
                    "effStat_idipMinus" : ".*effStat.*_idip.*q0",
                    "effStat_recoPlus"  : ".*effStat.*_reco.*q1",
                    "effStat_recoMinus" : ".*effStat.*_reco.*q0",
                    "effStat_trackingPlus"  : ".*effStat.*_tracking.*q1",
                    "effStat_trackingMinus" : ".*effStat.*_tracking.*q0",
                    "effSyst" : ".*effSyst.*",
                    "pdfAndAlphaS" : ".*pdf(\d+|.*AlphaS)",
}
if useSmoothSF:
    diffNuisanceDict["effStat_isoEffData"] = ".*effStat.*_iso_effData"
    diffNuisanceDict["effStat_isoEffMC"] = ".*effStat.*_iso_effMC"
else:
    diffNuisanceDict["effStat_isoSF"] = ".*effStat.*_iso"

    
print()

printText("PREFIT HISTOGRAMS")
command = f"python w-mass-13TeV/plotPrefitTemplatesWRemnants.py {combineInputFile} {mainOutdir}/plotPrefitTemplatesWRemnants/"
if isWlike:
    command += " --wlike"
if not skipHistograms:
    print()
    safeSystem(command, dryRun=dryRun)


printText("HISTOGRAM SYST RATIOS")
# the command below does charge plus, but it is just for example
processes = "Zmumu" if isWlike else "Wmunu,Fake"
command  = f"python w-mass-13TeV/makeSystRatios.py {combineInputFile} {mainOutdir}/makeSystRatios/"
command += f" -p {processes} -c plus --plotStat"
if not skipSystRatios:
    for key,value in systRatiosDict.items():
        trueCommand = command + f" --systPostfix {key} -s '{value}'"
        print()
        safeSystem(trueCommand, dryRun=dryRun)
                                                        
    
## Now produce plots from fit results for data or Asimov (toys to be implemented)
for fit in fits:

    if skipData and fit == "Data": continue
    if onlyData and fit != "Data": continue
    
    typedir = "hessian" if fit == "Asimov" else "data"

    printText(f"Running plots for {fit} fit")
    fitresultsFile = f"{mainInputPath}/{subFolder}/fit/{typedir}/fitresults_123456789_{fit}_bbb1_cxs0.root"

    ########################################
    printText("IMPACTS")
    statImpact = "0.0710" if isWlike else "0.0230"
    # TODO: add cases for single charge (change input to --set-stat and use "--postfix charge" not to overwrite plots)
    command = f"python w-mass-13TeV/makeImpactsOnMW.py {fitresultsFile} -o {mainOutdir}/makeImpactsOnMW/"
    command += f" --set-stat {statImpact} --showTotal --scaleToMeV"
    if not skipImpacts:
        print()
        safeSystem(command, dryRun=dryRun)
        
    ########################################
    printText("NUISANCES AND COSTRAINTS")
    yAxisSetting = " --y-setting -1.0 -0.5 0 0.5 1.0" if fit == "Asimov" else " --y-setting -5.0 -3.0 0 3.0 5.0"
    command  = f"python w-mass-13TeV/diffNuisances.py {fitresultsFile} -o {mainOutdir}/diffNuisances/ "
    command += f" -a --format html --type hessian  --suffix {fit} {yAxisSetting}"
    if not skipNuisances:
        for key,value in diffNuisanceDict.items():
            trueCommand = command + f" --uniqueString {key} -p '{value}'"
            print()
            safeSystem(trueCommand, dryRun=dryRun)

    ########################################
    printText("PREFIT AND POSTFIT HISTOGRAMS (YIELDS AND UNCERTAINTIES)")
    command  = f"python w-mass-13TeV/postFitHistograms.py {fitresultsFile} -o {mainOutdir}/postFitHistograms/ "
    command += f" --suffix postVFP -l 16.8 --no2Dplot" # remove --no2Dplot to add all 2D histograms
    if not skipPostfitHistograms:
        print()
        safeSystem(command, dryRun=dryRun)
        
print()
print("THE END!")
print()
