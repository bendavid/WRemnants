import os
import glob
import pdb
import argparse
import json
import uproot 
import numpy as np
import pandas as pd

os.sys.path.append(os.path.expandvars('/home/d/dwalter/WRemnants/'))

from utilities import logging
logger = logging.setup_logger(__file__, 3, True)

# Write LATEX table for results of biastests.
#   The shift in the mass parameter is the bias
# This script is step 3 and done within WRemnants singularity, it collects the biastests results and creates a latex table

parser = argparse.ArgumentParser()
parser.add_argument("-i","--input", type=str, help="json file with paths to fit results")
args = parser.parse_args()

outDir = "/".join(args.input.split("/")[:-1])

with open(args.input,"r") as rfile:
    results = json.load(rfile)

def readImpacts(rtfile, group, sort=True, add_total=True, stat=0.0):
    histname = "nuisance_group_impact_nois" if group else "nuisance_impact_nois"
    impacts = rtfile[histname].to_hist()
    labels = np.array([impacts.axes[1].value(i) for i in range(impacts.axes[1].size)])
    total = rtfile["fitresults"][impacts.axes[0].value(0)+"_err"].array()[0]
    impacts = impacts.values()[0,:]
    if sort:
        order = np.argsort(impacts)
        impacts = impacts[order]
        labels = labels[order]
    if add_total:
        impacts = np.append(impacts, total)
        labels = np.append(labels, "Total")

    if stat > 0:
        idx = np.argwhere(labels == "stat")
        impacts[idx] = stat

    return impacts,labels

def read_result(rootfile):
    logger.info(f"read {rootfile}")
    with uproot.open(rootfile) as rtfile:

        impacts_group = {l: i for i, l in zip(*readImpacts(rtfile, True))}

        pull_mass = rtfile["fitresults"]["massShift100MeV"].array()[0]

        uncertainty_pdf = impacts_group[f"pdf{nominal.upper()}"]
        uncertainty_total = impacts_group["Total"]

        return pull_mass, uncertainty_pdf, uncertainty_total

nominals = []
pseudos = []
channels = []
masses = []
uncertainties_pdf = []
uncertainties_tot = []

for nominal, r_n in results.items():

    for pseudodata, r_np in r_n.items():

        for channel, filename in r_np.items():

            m, updf, ut = read_result(filename)

            nominals.append(nominal.upper())
            pseudos.append(pseudodata.upper())
            channels.append(channel)
            masses.append(m)
            uncertainties_pdf.append(updf)
            uncertainties_tot.append(ut)

data = pd.DataFrame({"nominal":nominals, "pseudo":pseudos, "channel":channels, "mass":masses, "unc_pdf": uncertainties_pdf, "unc_tot":uncertainties_tot})

for channel, df in data.groupby("channel"):

    df.sort_values(by=["nominal","pseudo"])

    pseudo = list(sorted(set(df["pseudo"].values)))

    outfile=f"{outDir}/result_table_{channel}.txt"
    logger.info(f"write {outfile}")
    with open(outfile, "w") as outfile:

        outfile.write(r"\begin{table}" +"\n")
        outfile.write(r"\topcaption{\label{table:pullsPDFUnc}"+"\n")
        outfile.write(r"""  Pulls table for fitting $(p_T^{\ell},\eta^{\ell})$ in the """+channel+r""" channel with only pdf and bin-by-bin statistical uncertainties. 
  Entries read (pull on the $m_W$ central value) $\pm$ (pdf uncertainty on $m_W$).}"""+"\n")
        outfile.write(r"\centering"+"\n")

        columns = "l|"
        columns += "".join(["c" for c in range(len(pseudo))])
        outfile.write(r"\begin{tabular}{"+columns+"}"+"\n")
        
        outfile.write("  Model + & \multicolumn{"+str(len(pseudo))+"}{c}{Pseudodata} " + r" \\"+"\n")
        outfile.write("  Uncertainty & " + " & ".join(pseudo) + r" \\"+"\n")

        outfile.write(r"  \hline "+"\n")

        for nominal, df_n in df.groupby("nominal"):
            entries = []
            for p in pseudo:
                m, u = df_n.loc[df_n["pseudo"] == p][["mass","unc_pdf"]].values[0]
                m = round(100*m,1)
                u = round(100*u,1)
                colorstring = "\cellcolor{red!25}" if abs(m) > u else ""    # highlight background color of cell if uncertainty does not cover mass shift
                entries.append(f"{colorstring} ${m} \pm {u}$")

            outfile.write(f"  {nominal} & " + " & ".join(entries) + r" \\"+"\n")


        outfile.write(r"  \end{tabular}"+"\n")
        outfile.write(r"\end{table} "+"\n")
