import os
import glob
import pdb
import json
os.sys.path.append(os.path.expandvars('/home/d/dwalter/WRemnants/'))

from utilities import logging
logger = logging.setup_logger(__file__, 3, True)

# Perform bias tests for different pseudodata.
#   The shift in the mass parameter is the bias
# This script is step 2 and done within the cmssw-cc7 singularity and performs the fits 

parser = argparse.ArgumentParser()
parser.add_argument("-i","--input", type=str, help="input path to combine subfolders")
args = parser.parse_args()

histmaker = "mw_with_mu_eta_pt"
fittype = "WMass_pt_eta"   

combineDir = args.input

def EXE(command):
    logger.info(command) 
    os.system(command)  # for testing comment out this line

def get_card(name):
    card = glob.glob(name)
    if len(card) == 0:
        logger.error(f"Card {name} does not exist")
        return None
    elif len(card) > 1:
        logger.error(f"Multiple cards found for {name}")
        return None
    else:
        return card[0]


for subdir in glob.glob(f"{combineDir}/*_vs_*/{fittype}/"):

    nominal, pseudodata = subdir.split("/")[-2].split("_vs_")
    nominal = nominal.split("_")[-1]

    if nominal == pseudodata:
        continue

    if nominal not in results.keys():
        results[nominal] = {}
    if pseudodata not in results[nominal].keys():
        results[nominal][pseudodata] = {}

    logger.info(f"Now at {subdir}")
    
    card_minus = get_card(f"{subdir}/*_minus.txt")
    card_plus = get_card(f"{subdir}/*_plus.txt")
    card_combined = card_minus.replace("_minus.txt",".txt")

    if not os.path.isfile(card_combined):
        # merge cards
        # first move them into the current directory, otherwise the merging does not work, then move them back
        EXE(f"mv {card_minus} {card_plus} ./")
        EXE("combineCards.py "+card_minus.split("/")[-1]+" "+card_plus.split("/")[-1]+" > "+card_combined.split("/")[-1])
        EXE("mv "+card_combined.split("/")[-1]+" "+card_combined)
        EXE("mv "+card_minus.split("/")[-1]+" "+card_minus)
        EXE("mv "+card_plus.split("/")[-1]+" "+card_plus)
    
    if nominal not in results[nominal].keys():
        results[nominal][nominal] = {}
    
    results[nominal][pseudodata] = {}

    for card, channel in (
        (card_minus, "minus"), 
        (card_plus, "plus"),
        (card_combined, "combined")
        ):
        logger.info(f"Now at {channel}") 

        input_hdf5 = card.replace(".txt",".hdf5")

        if not os.path.isfile(input_hdf5):
            # convert .txt into .hdf5
            EXE(f"text2hdf5.py {card} --X-allow-no-signal")

        fitresult = card.replace(".txt","_fitresult.root")

        if not os.path.isfile(fitresult):
            # run combine fit to pseudodata
            EXE(f"combinetf.py {input_hdf5} --doImpacts --binByBinStat -o {fitresult}")
        else:
            logger.info("Existing fitresult file found")

        # make asimov fit if it was not done already
        if channel not in results[nominal][nominal].keys():
            fitresult_asimov = fitresult.replace(".root","_asimov.root")

            if not os.path.isfile(fitresult_asimov):
                # run combine fit to asimov data
                EXE(f"combinetf.py -t -1 {input_hdf5} --doImpacts --binByBinStat -o {fitresult_asimov}")
            else:
                logger.info("Existing fitresult file found")

            results[nominal][nominal][channel] = fitresult_asimov

        results[nominal][pseudodata][channel] = fitresult

with open(f"{combineDir}/results.json","w") as rfile:
    json.dump(results, rfile, indent=4)
