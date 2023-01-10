# WRemnants

## Instructions:

Activate the singularity image (to be done every time before running code)
```bash
singularity run /cvmfs/unpacked.cern.ch/gitlab-registry.cern.ch/bendavid/cmswmassdocker/wmassdevrolling\:latest
```
    
Get the code (after forking from the central WMass repository)
```bash
MY_GIT_USER=$(git config user.github)
git clone --recurse-submodules git@github.com:$MY_GIT_USER/WRemnants.git
cd WRemnants/
git remote add upstream git@github.com:WMass/WRemnants.git
```

Get updates from the central repository (and main branch)
```bash
git pull --recurse-submodules upstream main
git push origin main
```
    
### Run the code
Source the setup script.
It will create some environment variables to ease access to some output paths:
 * WREM_BASE: it points to ./WRemnants/ where all the code is
 * COMBINE_STUDIES: main folder to store datacards and root files for combine
 * PLOTS: where to store analysis plots
These are used to define default output paths, but users can specify their favorite one for each script.
```bash
source WRemnants/setup.sh
```

Make histograms for WMass (for Wlike the script is __mz_wlike_with_mu_eta_pt.py__).
```
python WRemnants/scripts/histmakers/mw_with_mu_eta_pt.py --theory_corr scetlib --muonCorr none --no_recoil -p smoothSF_muonCorrNone_noRecoil
```
More options are loaded from **WRemnants/utilities/common.py**

Make the datacards for single charges and prepare the TH2 histograms for combinetf.
```
python WRemnants/scripts/combine/setupCombineWMass.py -i mw_with_mu_eta_pt_scetlibCorr_nnpdf31_smoothSF_muonCorrNone_noRecoil.pkl.lz4 -o smoothSF/muonCorr_none/scetlibCorr_nnpdf31/byHelicityPtCharge/ --skipOtherChargeSyst
```
The input file is the output of the previous step.
Add __--wlike__ if running a Wlike analysis.
The path specified with __-o__ is appended to the argument of option __-d__, whose default is \$COMBINE_STUDIES/WMass/ or \$COMBINE_STUDIES/ZMassWLike/ depending on the chosen analysis. One can also use __-d__ to pass any custom path (__-o__ is optional in that case)
 
Combine the datacards for single charges and run the fit (Asimov only)
```
python WRemnants/scripts/combine/fitManager.py -i $COMBINE_STUDIES/WMass/smoothSF/muonCorr_none/scetlibCorr_nnpdf31/byHelicityPtCharge/ --comb --skip-fit-data
```
Run the fits for single charges (Asimov only). These can be produced in the same output folder as the combination, since a postfix is automatically appended to the output card and fit result files.
```
python WRemnants/scripts/combine/fitManager.py -i $COMBINE_STUDIES/WMass/smoothSF/muonCorr_none/scetlibCorr_nnpdf31/byHelicityPtCharge/ --fit-single-charge --skip-fit-data [-c plus|minus|both]
```

**NOTE**: to run __fitManager.py__ one has to set a Centos 7 environment with __cmsssw-cc7__. Then, one has to activate __cmsenv__ from the folder where combine is installed (once the environment is set one can keep working from inside WRemnants).

Make plots on prefit histograms (output root file from __setupCombineWMass.py__) or the fit results (output of __fitManager.py__)
```
cd $WREM_BASE/scripts/analysisTools/
python runAllPlots.py
```
**NOTE**: Before running __runAllPlots.py__ make sure the settings inside it are adequate, especially the input/output folders.