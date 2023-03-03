import os
import glob
import pdb

os.sys.path.append(os.path.expandvars('/home/d/dwalter/WRemnants/'))

from utilities import logging
logger = logging.setup_logger(__file__, 3, True)

# Perform bias tests on different pdf sets. Use one PDF set + uncertainties as nominal and fit the other pdf set as pseudodata.
#   The shift in the mass parameter is the bias
# This script is step 3 and done within WRemnants singularity, it collects the biastests results and creates a latex table

pdfsets = ["nnpdf31", "nnpdf40", "msht20", "pdf4lhc21", "ct18"]

histmaker = "mw_with_mu_eta_pt"
fittype = "WMass_pt_eta"   

combineDir = "/scratch/dwalter/CombineStudies/230302_studies_pdfs/"

webdir = "/eos/home-d/dwalter/www/WMassAnalysis/230302_studies_pdfs/"

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

# def readImpacts(rtfile, group, sort=True, add_total=True, stat=0.0):
#     histname = "nuisance_group_impact_nois" if group else "nuisance_impact_nois"
#     impacts = rtfile[histname].to_hist()
#     labels = np.array([impacts.axes[1].value(i) for i in range(impacts.axes[1].size)])
#     total = rtfile["fitresults"][impacts.axes[0].value(0)+"_err"].array()[0]
#     impacts = impacts.values()[0,:]
#     if sort:
#         order = np.argsort(impacts)
#         impacts = impacts[order]
#         labels = labels[order]
#     if add_total:
#         impacts = np.append(impacts, total)
#         labels = np.append(labels, "Total")

#     if stat > 0:
#         idx = np.argwhere(labels == "stat")
#         impacts[idx] = stat

#     return impacts,labels

# def read_result(rootfile):
#     with uproot.open(args.inputFile) as rtfile:

#         impacts = {l: i for i, l in zip(readImpacts(rtfile, False))}
#         impacts_group = {l: i for i, l in zip(readImpacts(rtfile, True))}

#         pull_mass = getattr(rtfile[treename], "MassShift100MeV")

#         uncertainty = impacts_group[f"pdf{nominal.upper()}"]
#         # uncertainty = impacts_group[f"Total"]

#         return pull_mass, uncertainty

results = {}
for nominal in pdfsets:
    results[nominal] = {}

    for pseudodata in pdfsets:
        if nominal == pseudodata:
            asimov=True

        for subdir in glob.glob(f"{combineDir}/*_{nominal}_vs_{pseudodata}/{fittype}/"):
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
                    pdb.set_trace()
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

                    # read out results
                    # m, u = read_result(fitresult_asimov)
                    #     "mass": m,
                    #     "uncertainty": u
                    # }

                # read out fit results
                # m, u = read_result(fitresult)
                results[nominal][pseudodata][channel] = fitresult
                #     "mass": m,
                #     "uncertainty": u
                # }

            exit()

                # impactfile = f"{webdir}/impacts_{nominal}_vs_{pseudodata}_{channel}.html"
                # if not os.path.isfile(impactfile):
                #     EXE(f"python3 scripts/combine/pullsAndImpacts.py -f {fitresult} -s impact output -o {impactfile}")

                # impactfile_asimov = f"{webdir}/impacts_{nominal}_{channel}.html"
                # if not os.path.isfile(impactfile_asimov):
                #     EXE(f"python3 scripts/combine/pullsAndImpacts.py -f {fitresult} -s impact output -o {impactfile_asimov}")


