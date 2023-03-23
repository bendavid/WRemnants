import os
import glob
import pdb
import argparse
import json
import uproot 
import numpy as np
import pandas as pd

os.sys.path.append(os.path.expandvars('/home/d/dwalter/WRemnants/'))

from utilities import logging, input_tools
logger = logging.setup_logger(__file__, 3, True)

# Write LATEX table for results of biastests.
#   The shift in the mass parameter is the bias
# This script is step 3 and done within WRemnants singularity, it collects the biastests results and creates a latex table

parser = argparse.ArgumentParser()
parser.add_argument("-i","--input", type=str, help="json file with paths to fit results")
parser.add_argument("-m","--mode", type=str, default="tot", help="Uncertainty to consider", choices=["pdf","scale"])
parser.add_argument("-f","--folder", type=str, default="Impacts", help="output folder for impact plots")
args = parser.parse_args()

outDir = "/".join(args.input.split("/")[:-1])

do_impacts = False
group = True
group_str = "-g" if group else ""

output_dir_impacts = f"/eos/home-d/dwalter/www/WMassAnalysis/{args.folder}"
if do_impacts and not os.path.isdir(output_dir_impacts):
    os.mkdir(output_dir_impacts)

def EXE(command):
    logger.info(command) 
    os.system(command)  # for testing comment out this line

def read_result(rootfile, nominal):
    logger.info(f"read {rootfile}")
    with uproot.open(rootfile) as rtfile:

        impacts_group = {l: i for i, l in zip(*input_tools.readImpacts(rtfile, True))}

        pull_mass = rtfile["fitresults"]["massShift100MeV"].array()[0]

        if args.mode == "scale":
            uncertainty = impacts_group["Total"]
        elif args.mode == "pdf":
            uncertainty = impacts_group[f"pdf{nominal.upper()}"]
        

        return pull_mass, uncertainty


with open(args.input,"r") as rfile:
    results = json.load(rfile)

nominals = []
pseudos = []
channels = []
masses = []
uncertainties = []

update_result = False
for nominal, r_n in results.items():

    for pseudodata, r_np in r_n.items():

        for channel, r_c in r_np.items():

            if isinstance(r_c, dict) and "mass" in r_c.keys() and "unc" in r_c.keys():
                m = results[nominal][pseudodata][channel]["mass"]
                u = results[nominal][pseudodata][channel]["unc"]
                filename = results[nominal][pseudodata][channel]["filename"]
            elif isinstance(r_c, str):
                filename = r_c
                m, u = read_result(filename, nominal)

                results[nominal][pseudodata][channel] = {
                    "filename" : filename,
                    "mass": m,
                    "unc":u
                }
                update_result = True
            
            
            output_impacts = f"{output_dir_impacts}/{channel}_{nominal}_vs_{pseudodata}.html"

            if do_impacts and not os.path.isfile(output_impacts):
                
                EXE(f"python3 scripts/combine/pullsAndImpacts.py -f {filename} {group_str} -s impact output -o {output_impacts}")

            if args.mode == "pdf":
                nominals.append(nominal.upper())
                pseudos.append(pseudodata.upper())
            else:
                nominals.append(nominal)
                pseudos.append(pseudodata)

            channels.append(channel)
            masses.append(m)
            uncertainties.append(u)



if update_result:
    # save result.json so that next time the information does not have to be read again
    with open(args.input,"w") as rfile:
        json.dump(results, rfile, indent=4)


### Make latex table

# some settings
boson_str = "W"
fit_str = "\P"+boson_str+" $(\pt^{\ell},\eta^{\ell})$"

if args.mode == "scale":
    mode_str = "muon momentum scale and resoltion"
    uncertainty_str1 = "with the full uncertainty model."
    uncertainty_str2 = " uncertainty on $m_\P"+boson_str+"$"

elif args.mode == "pdf":
    mode_str = "pdf sets"
    uncertainty_str1 = "only pdf and bin-by-bin statistical uncertainties."
    uncertainty_str2 = "pdf uncertainty on $m_\P"+boson_str+"$"


data = pd.DataFrame({"nominal":nominals, "pseudo":pseudos, "channel":channels, "mass":masses, "unc": uncertainties})

if args.mode == "pdf":
    for channel, df in data.groupby("channel"):

        df.sort_values(by=["nominal","pseudo"])

        pseudo = list(sorted(set(df["pseudo"].values)))
        nominal = list(sorted(set(df["nominal"].values)))

        pseudo = list(filter(lambda x: x in pseudo, nominal)) + list(filter(lambda x: x not in nominal, pseudo))

        outfile=f"{outDir}/result_table_{channel}.txt"
        logger.info(f"write {outfile}")
        with open(outfile, "w") as outfile:

            outfile.write(r"\begin{table}" +"\n")
            outfile.write(r"\topcaption{\label{table:pulls_pdf_"+channel+"}"+"\n")
            outfile.write(r""" Pulls table for different pdf sets. The fit is performed on """+fit_str+""" in the """+channel+r""" channel with """+uncertainty_str1+r""" 
    Entries read (pull on the $m_\P"""+boson_str+"""$ central value) $\pm$ (""" +uncertainty_str2+ r""").}"""+"\n")
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
                    df_p = df_n.loc[df_n["pseudo"] == p][["mass", f"unc" ]]        
                    if len(df_p) == 1:
                        m, u = df_p.values[0]
                        m = round(100*m,1)
                        u = round(100*u,1)
                        colorstring = "\cellcolor{red!25}" if abs(m) > u else ""    # highlight background color of cell if uncertainty does not cover mass shift
                        entries.append(f"{colorstring} ${m} \pm {u}$")
                    else:
                        entries.append(r" \NA ")

                outfile.write(f"  {nominal} & " + " & ".join(entries) + r" \\"+"\n")


            outfile.write(r"  \end{tabular}"+"\n")
            outfile.write(r"\end{table} "+"\n")

if args.mode == "scale":

    translate = {
        "smearing": "smeared",
        "biasCalibrationA": "bias A",
        "biasCalibrationM": "bias M",
        "biasCalibrationParameterized": "bias A,M",
        "biasCalibrationBinned": "bias $\pt,\eta$",
        "smearing_biasCalibrationA": "+ bias A",
        "smearing_biasCalibrationM": "+ bias M",
        "smearing_biasCalibrationParameterized": "+ bias A,M",
        "smearing_biasCalibrationBinned": "+ bias $\pt,\eta$"
    }

    for nominal, df in data.groupby("nominal"):

        closure_str = "different muon momentum scale closures" if nominal=="smearing" else "muon momentum scale resolution"

        df.sort_values(by=["pseudo"])

        pseudo = list(sorted(set(df["pseudo"].values)))
        pseudo = list(filter(lambda x: x in pseudo, nominal)) + list(filter(lambda x: x not in nominal, pseudo))

        outfile=f"{outDir}/result_table_{nominal}.txt"
        logger.info(f"write {outfile}")
        with open(outfile, "w") as outfile:

            outfile.write(r"\begin{table}" +"\n")
            outfile.write(r"\topcaption{\label{table:pulls_scale_"+nominal+"}"+"\n")
            outfile.write(r""" Pulls table for """+closure_str+""". The fit is performed on """+fit_str+""" in different channels with """+uncertainty_str1+r""" 
    Entries read (pull on the $m_\P"""+boson_str+"""$ central value) $\pm$ (""" +uncertainty_str2+ r""").}"""+"\n")
            outfile.write(r"\centering"+"\n")

            columns = "l|"
            columns += "".join(["c" for c in range(len(pseudo))])
            outfile.write(r"\begin{tabular}{"+columns+"}"+"\n")
            
            outfile.write("  Model + & \multicolumn{"+str(len(pseudo))+"}{c}{Pseudodata} " + r" \\"+"\n")
            outfile.write("  Uncertainty & " + " & ".join([translate.get(p,p) for p in pseudo]) + r" \\"+"\n")

            outfile.write(r"  \hline "+"\n")

            for channel, df_n in df.groupby("channel"):
                entries = []
                for p in pseudo:         
                    df_p = df_n.loc[df_n["pseudo"] == p][["mass", f"unc" ]]        
                    if len(df_p) == 1:
                        m, u = df_p.values[0]
                        m = round(100*m,1)
                        u = round(100*u,1)
                        colorstring = "\cellcolor{red!25}" if abs(m) > u else ""    # highlight background color of cell if uncertainty does not cover mass shift
                        entries.append(f"{colorstring} ${m} \pm {u}$")
                    else:
                        entries.append(r" \NA ")

                outfile.write(f"  {translate.get(nominal,nominal)} ({channel}) & " + " & ".join(entries) + r" \\"+"\n")


            outfile.write(r"  \end{tabular}"+"\n")
            outfile.write(r"\end{table} "+"\n")