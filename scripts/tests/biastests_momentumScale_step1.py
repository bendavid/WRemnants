import os
from utilities import common

logger = common.setup_logger(__file__, 3, True)

# Perform bias tests on different momentum scale and resolution corrections and closure files. 
#   Use one setting + uncertainties as nominal and fit the other setting as pseudodata.
#   The shift in the mass parameter is the bias
# This script is step 1 and done within the WRemnants singularity and produces the histograms and combine inputs
#   A second script, `biastests_momentumScale_step2.py` , has to be run in the cmsssw-7 environment and the fits are performed and the information is read out

histmaker = "mw_with_mu_eta_pt"

histDir = "/scratch/dwalter/results_histmaker/230303_biastests_momentumScale/"
combineDir = "/scratch/dwalter/CombineStudies/230303_biastests_momentumScale/"

nTreads = 128

# common options (apart from --theory_corr)
theory_corr = "scetlib"
options = "--no-recoil"

if not os.path.isdir(histDir):
    os.mkdir(histDir)

if not os.path.isdir(combineDir):
    os.mkdir(combineDir)

def EXE(command):
    logging.info(command) 
    os.system(command)  # for testing comment out this line

def make_appendix(name):
        
    name = "{options} {name}"
    parts = []
    for p in [p.replace("-"," ") for p in name.split("--")] :
        ps = p.split(" ")
        ps = "".join([ps[0],]+[p.capitalize() for p in ps[1:]])
        parts.append(ps)

    if theory_corr != "":
        theory_append = theory_corr+"Corr_"

    return theory_append + "_".join(parts)


for nominal, pseudodata in (
    ("", "--bias-calibration binned"),
    ("", "--bias-calibration parameterized"),
    ("", "--smearing"),
    ):
    if nominal == pseudodata:
        logging.warning("Now at nominal and pseudodata are the same, this would be just an asimov fit! continue with next setup")
        continue

    logging.info("Now at {nominal} vs {pseudodata}")

    append_nominal = make_appendix(nominal)
    append_pseudodata = make_appendix(pseudodata)

    file_nominal = f"{histDir}/{histmaker}{append_nominal}.hdf5"
    file_pseudodata = f"{histDir}/{histmaker}{append_pseudodata}.hdf5"

    if not os.path.isfile(file_nominal):
        # run histmaker for nominal
        EXE(f"python3 scripts/histmakers/{histmaker}.py -j {nTreads} -o {histDir} --theory_corr {theory_corr} {options} -p {append_nominal}")
    else:
        logging.info(f"Found file for nominal {nominal}")

    if not os.path.isfile(file_pseudodata):
        # run histmaker for pseudodata
        EXE(f"python3 scripts/histmakers/{histmaker}.py -j {nTreads} -o {histDir} --theory_corr {theory_corr} {options} -p {append_pseudodata}")
    else:
        logging.info(f"Found file for pseudodata {pseudodata}")

    # make combine input
    dir_combine = f"{combineDir}/{append_nominal}_vs_{append_pseudodata}"

    if os.path.isdir(dir_combine):
        logging.warning("The combine file for {dir_combine} already exists, continue with the next one!")
        continue

    EXE(f"python3 scripts/combine/setupCombineWMass.py -o {dir_combine} -i {file_nominal} --pseudodata-file {file_pseudodata} --pseudoData nominal")

    exit()

