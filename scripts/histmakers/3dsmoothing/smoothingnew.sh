#!/bin/bash

#SCRIPT TO HELP WITH THE PRODUCTION OF THE HISTOGRAMS TO DEFINE THE 3D EFFICIENCIES INTEGRATED OVER THE W UT DISTRIBUTION. YOU CAN OF COURSE CHANGE theory_corr, no_recoil etc...

NAME="20bins"

python WRemnants/scripts/histmakers/mw_with_mu_eta_pt.py --theory_corr scetlib --muonCorr none --no_recoil -p test${NAME} --binnedScaleFactors --vqtTest --vqtTestReal --vqtTestStep 2 --vqtTestCorrectionStep 2
