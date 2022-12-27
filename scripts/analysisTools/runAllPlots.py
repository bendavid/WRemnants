#!/bin/env pyth

### Wrapper to run all plotting scripts in one go out of the fit results
### This runs the minimal command, which is usually sufficient
### Each script has additional options for more customization, please check them

import ROOT, os, sys, re, array
sys.path.append(os.getcwd() + "/plotUtils/")
from plotUtils.utility import safeSystem

# general settings
dryRun = 1
isWlike = 0
skipData = 1 # or set fits = ["Asimov"]
onlyData = 0 # or set fits = ["Data"]
fits = ["Asimov", "Data"]

# what to plot
skipHistograms = 0 # prefit histograms, doesn't require having run the fit
skipImpacts = 0
skipNuisances = 0
skipSystRatios = 0
skipPostfitHistograms = 0 # prefit and postfit histograms, from fitresults.root

# input and output folders
# TODO: need a general way so that every user can use the same logic for paths
mainInputPath = "/scratch/mciprian/CombineStudies/WMass/scetlibCorr_nnpdf31_testSF_trashdebug/byHelicityPtCharge_correlateEffStatIsoByCharge/" # contains the TH2 histograms for combine
subFolder = "nominal" # contains the final cards and fit results
mainOutdir = "/eos/user/m/mciprian/www/WMassAnalysis/fromMyWremnants/Wmass_fit/TEST_PLOT_SCRIPTS/" # where to store plots

## These should not be touched
what = "WLike" if isWlike else "WMass"
combineInputFile = f"{mainInputPath}/{what}CombineInput.root" 


##############################
# utility functions used here
def printText(text):
    print("")
    print("="*30)
    print(text)
    print("")


##############################
# to customize some specific plots
#
#
# postfix for plot name and regular expression to pick histogram variations to plot
# there are infinite combinations, so this is mainly an example for Up variations
systRatiosDict = {"effStatTnP_trigger_eta20_plus"  : ".*effStatTnP.*_trigger_eta20.*q1.*Up",
                  "pdfsAndAlphaS"  : ".*pdf(12|25|50|.*AlphaS).*Up",
}

#
# unique string for output plot and regular expression for nuisances used to plot pulls and constraints
diffNuisanceDict = {"effStat_triggerPlus"  : ".*effStatTnP.*_trigger.*q1",
                    "effStat_triggerMinus" : ".*effStatTnP.*_trigger.*q0",
                    "effStat_isoEffData"   : ".*effStatTnP.*_iso_effData",
                    "effSyst" : ".*effSystTnP_",
                    "pdfAndAlphaS" : ".*pdf(\d+|.*AlphaS)",
}

    
print()

printText("PREFIT HISTOGRAMS")
command = f"python w-mass-13TeV/plotPrefitTemplatesWRemnants.py {combineInputFile} {mainOutdir}/plotPrefitTemplatesWRemnants/"
if isWlike:
    command += " --wlike"
if not skipHistograms:
    print()
    safeSystem(command, dryRun=dryRun)


printText("HISTOGRAM SYST RATIOS")
# by default it does charge plus, it is just for example
processes = "Zmumu" if isWlike else "Wmunu,Fake"
command  = f"python w-mass-13TeV/makeSystRatios.py {combineInputFile} {mainOutdir}/makeSystRatios/"
command += f" -p {processes} -c plus --plotStat"
if not skipSystRatios:
    for key,value in systRatiosDict.items():
        trueCommand = command + f" --systPostfix {key} -s '{value}'"
        print()
        safeSystem(trueCommand, dryRun=dryRun)
                                                        
    
## Now plots from fit results
for fit in fits:

    if skipData and fit == "Data": continue
    if onlyData and fit != "Data": continue
    
    typedir = "hessian" if fit == "Asimov" else "data"

    printText(f"Running plots for {fit} fit")

    fitresultsFile = f"{mainInputPath}/{subFolder}/fit/{typedir}/fitresults_123456789_{fit}_bbb1_cxs0.root"

    ##
    printText("IMPACTS")
    # TODO: add cases for single charge (change input to --set-stat and use "--postfix charge" not to overwrite plots)
    command = f"python w-mass-13TeV/makeImpactsOnMW.py {fitresultsFile} -o {mainOutdir}/makeImpactsOnMW/ --set-stat 0.0230 --showTotal --scaleToMeV"
    if not skipImpacts:
        print()
        safeSystem(command, dryRun=dryRun)

    ##
    printText("NUISANCES AND COSTRAINTS")
    yAxisSetting = " --y-setting -1.0 -0.5 0 0.5 1.0" if fit == "Asimov" else " --y-setting -5.0 -3.0 0 3.0 5.0"
    command  = f"python w-mass-13TeV/diffNuisances.py {fitresultsFile} -o {mainOutdir}/diffNuisances/ "
    command += f" -a --format html --type hessian  --suffix {fit} {yAxisSetting}"
    if not skipNuisances:
        for key,value in diffNuisanceDict.items():
            trueCommand = command + f" --uniqueString {key} -p '{value}'"
            print()
            safeSystem(trueCommand, dryRun=dryRun)

    ##
    printText("PREFIT AND POSTFIT HISTOGRAMS (YIELDS AND UNCERTAINTIES)")
    command  = f"python w-mass-13TeV/postFitHistograms.py {fitresultsFile} -o {mainOutdir}/postFitHistograms/ "
    command += f" --suffix postVFP -l 16.8 --no2Dplot" # remove --no2Dplot to add all 2D histograms
    if not skipPostfitHistograms:
        print()
        safeSystem(command, dryRun=dryRun)
        
print()
print("THE END!")
print()
