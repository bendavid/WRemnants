#!/bin/bash

NAME="utdistrofull"

python WRemnants/scripts/histmakers/mw_with_mu_eta_pt.py --theory_corr scetlib --muonCorrMC none --muonCorrData none --no_recoil -p idip${NAME} --binnedScaleFactors --vqtTest --vqtTestReal --vqt3dsmoothing --vqtTestStep 0 --vqtTestCorrectionStep 0
python WRemnants/scripts/histmakers/mw_with_mu_eta_pt.py --theory_corr scetlib --muonCorrMC none --muonCorrData none --no_recoil -p iso${NAME} --binnedScaleFactors --vqtTest --vqtTestReal --vqt3dsmoothing --vqtTestStep 2 --vqtTestCorrectionStep 2
