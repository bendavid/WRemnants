import os
from utilities import logging
import argparse 

logger = logging.setup_logger(__file__, 3, True)

# Perform bias tests on different (pseudo) data sets. Use one dataset + uncertainties as nominal and fit the other set as pseudodata.
#   The shift in the mass parameter is the bias
# This script is step 1 and done within the WRemnants singularity and produces the histograms and combine inputs
#   A second script, `biastests_step2.py` , has to be run in the cmsssw-7 environment and the fits are performed and the information is read out

parser = argparse.ArgumentParser()
parser.add_argument("-m","--mode", type=str, help="what kind of biastests", choices=["pdf","scale"])
args = parser.parse_args()


if args.mode == "pdf":
    argument = "--pdfs"

    nominals = ["nnpdf31", "nnpdf40", "msht20", "pdf4lhc21", "ct18"]
    pseudos = ["nnpdf31", "nnpdf40", "msht20", "pdf4lhc21", "ct18"]

    # freeze_uncertainties = ("reducedUncertainties", "-x '.*' -k 'pdf*|alphaS*|mass*'") # freeze all nuisances except pdf and alphaS
    freeze_uncertainties = ("fullUncertainties", "")

    histDir = "/scratch/dwalter/results_histmaker/230305_studies_pdfs/"
    combineDir = f"/scratch/dwalter/CombineStudies/230305_studies_pdfs_{freeze_uncertainties[0]}/"

elif args.mode == "scale":
    argument = ""

    nominals = ["default",]
    pseudos = ["--bias-calibration binned", "--bias-calibration parameterized", "--smearing"]

    freeze_uncertainties = ("fullUncertainties", "")

    histDir = "/scratch/dwalter/results_histmaker/230305_biastests_momentumScale/"
    combineDir = f"/scratch/dwalter/CombineStudies/230305_biastests_momentumScale_{freeze_uncertainties[0]}/"


histmaker = "mw_with_mu_eta_pt"

nTreads = 128

options = [
    "--theory_corr scetlib",
    "--no_recoil"
]
options = " ".join(options)

if not os.path.isdir(histDir):
    os.mkdir(histDir)

if not os.path.isdir(combineDir):
    os.mkdir(combineDir)

def EXE(command):
    logger.info(command) 
    os.system(command)  # for testing comment out this line

def make_appendix(name):
        
    parts = []
    for p in [pa.replace("-"," ").replace("_"," ") for pa in name.split("--")]:
        if p=="":
            continue
        ps = p.split(" ")
        ps = "".join([ps[0],]+[p.capitalize() for p in ps[1:]])
        parts.append(ps)
    
    return "_".join(parts)

for nominal in nominals:
    
    for pseudodata in pseudos:
        if nominal == pseudodata:
            continue

        logger.info(f"Now at {nominal} vs {pseudodata}")
        
        arg_nominal = nominal if nominal != "default" else ""
        str_nominal = make_appendix(nominal)

        arg_pseudodata = pseudodata if pseudodata != "default" else ""
        str_pseudodata = make_appendix(pseudodata)

        file_nominal = f"{histDir}/{histmaker}_scetlibCorr_noRecoil_{str_nominal}.hdf5"
        file_pseudodata = f"{histDir}/{histmaker}_scetlibCorr_noRecoil_{str_pseudodata}.hdf5"

        if not os.path.isfile(file_nominal):
            # run histmaker for nominal
            EXE(f"python3 scripts/histmakers/{histmaker}.py -j {nTreads} -o {histDir} {options} {argument} {arg_nominal} -p noRecoil_{str_nominal}")
        else:
            logger.info(f"Found file for nominal {nominal}")

        if not os.path.isfile(file_pseudodata):
            # run histmaker for pseudodata
            EXE(f"python3 scripts/histmakers/{histmaker}.py -j {nTreads} -o {histDir} {options} {argument} {arg_pseudodata} -p noRecoil_{str_pseudodata}")
        else:
            logger.info(f"Found file for pseudodata {pseudodata}")

        # make combine input        
        dir_combine = f"{combineDir}/{str_nominal}_vs_{str_pseudodata}"

        if os.path.isdir(dir_combine):
            logger.warning(f"The combine file for {dir_combine} already exists, continue with the next one!")
            continue

        EXE(f"python3 scripts/combine/setupCombineWMass.py -o {dir_combine} -i {file_nominal} --pseudodata-file {file_pseudodata} --pseudoData nominal {freeze_uncertainties[1]}")
    

