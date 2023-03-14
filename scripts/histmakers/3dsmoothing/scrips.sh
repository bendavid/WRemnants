#!/bin/bash

#SCRIPT TO HELP WITH THE PRODUCTION OF THE HISTOGRAMS TO DEFINE THE 3D EFFICIENCIES INTEGRATED OVER THE W UT DISTRIBUTION. YOU CAN OF COURSE CHANGE theory_corr, no_recoil etc...

NAME="redo"

python WRemnants/scripts/histmakers/mw_with_mu_eta_pt.py --theory_corr scetlib --muonCorrMC none --muonCorrData none --no_recoil -p idip${NAME} --binnedScaleFactors --vqtTest --vqtTestReal --vqt3dsmoothing --vqtTestStep 0 --vqtTestCorrectionStep 0
python WRemnants/scripts/histmakers/mw_with_mu_eta_pt.py --theory_corr scetlib --muonCorrMC none --muonCorrData none --no_recoil -p trigger${NAME} --binnedScaleFactors --vqtTest --vqtTestReal --vqt3dsmoothing --vqtTestStep 1 --vqtTestCorrectionStep 1
python WRemnants/scripts/histmakers/mw_with_mu_eta_pt.py --theory_corr scetlib --muonCorrMC none --muonCorrData none --no_recoil -p iso${NAME} --binnedScaleFactors --vqtTest --vqtTestReal --vqt3dsmoothing --vqtTestStep 2 --vqtTestCorrectionStep 2
