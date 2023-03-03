import os
from utilities import logging

logger = logging.setup_logger(__file__, 3, True)

# Perform bias tests on different pdf sets. Use one PDF set + uncertainties as nominal and fit the other pdf set as pseudodata.
#   The shift in the mass parameter is the bias
# This script is step 1 and done within the WRemnants singularity and produces the histograms and combine inputs
#   A second script, `biastests_pdf_step2.py` , has to be run in the cmsssw-7 environment and the fits are performed and the information is read out

pdfsets = ["nnpdf31", "nnpdf40", "msht20", "pdf4lhc21", "ct18"]

histmaker = "mw_with_mu_eta_pt"

histDir = "/scratch/dwalter/results_histmaker/230302_studies_pdfs/"
combineDir = "/scratch/dwalter/CombineStudies/230302_studies_pdfs/"

nTreads = 128

options = [
    "--filterProcs 'data'", "--invert-filter", # skip data
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

for nominal in pdfsets:
    
    for pseudodata in pdfsets:
        if nominal == pseudodata:
            continue

        logger.info("Now at {nominal} vs {pseudodata}")
        
        file_nominal = f"{histDir}/{histmaker}_scetlibCorr_noRecoil_{nominal}.hdf5"
        file_pseudodata = f"{histDir}/{histmaker}_scetlibCorr_noRecoil_{pseudodata}.hdf5"

        if not os.path.isfile(file_nominal):
            # run histmaker for nominal
            EXE(f"python3 scripts/histmakers/{histmaker}.py -j {nTreads} -o {histDir} {options} --pdfs {nominal} -p noRecoil_{nominal}")
        else:
            logger.info(f"Found file for nominal {nominal}")

        if not os.path.isfile(file_pseudodata):
            # run histmaker for pseudodata
            EXE(f"python3 scripts/histmakers/{histmaker}.py -j {nTreads} -o {histDir} {options} --pdfs {pseudodata} -p noRecoil_{pseudodata}")
        else:
            logger.info(f"Found file for pseudodata {pseudodata}")

        # make combine input
        for freeze_string, freeze, in (
            ("reducedUncertainties", "-x '.*' -k 'pdf*|alphaS*|mass*'"),  # freeze all nuisances except pdf and alphaS
            ("fullUncertainties", "")
        ):     
            dir_combine = f"{combineDir}/{freeze_string}_{nominal}_vs_{pseudodata}"

            if os.path.isdir(dir_combine):
                logger.warning(f"The combine file for {dir_combine} already exists, continue with the next one!")
                continue

            EXE(f"python3 scripts/combine/setupCombineWMass.py -o {dir_combine} -i {file_nominal} --pseudodata-file {file_pseudodata} --pseudoData nominal {freeze}")
        
        exit()

