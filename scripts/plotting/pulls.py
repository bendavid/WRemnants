import mplhep as hep
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker

import uproot
import argparse
import os
import numpy as np
import pdb

from utilities import logging, output_tools
from wremnants import plot_tools

hep.style.use(hep.style.ROOT)

parser = argparse.ArgumentParser()
parser.add_argument("infile", help="Combine fitresult root file")
parser.add_argument("--asimov",  type=str, default=None, help="Optional combine fitresult root file from an asimov fit for comparison")
parser.add_argument("-o", "--outpath", type=str, default=os.path.expanduser("~/www/WMassAnalysis"), help="Base path for output")
parser.add_argument("-f", "--outfolder", type=str, default="./", help="Subfolder for output")
parser.add_argument("-p", "--postfix", type=str, help="Postfix for output file name")
parser.add_argument("--debug", action='store_true', help="Print debug output")
parser.add_argument("--groups", type=str, nargs="+", default=None, 
    choices=["scale", "prefire", "recoil", "pdf", "scetlib", "resum", "angularCoefficients", "eff_stat", "eff_syst"], 
    help="Define which plots to make")
parser.add_argument("--eoscp", action='store_true', help="Override use of xrdcp and use the mount instead")

args = parser.parse_args()

logger = logging.setup_logger("plotFitresult", 4 if args.debug else 3)

outdir = output_tools.make_plot_dir(args.outpath, args.outfolder)

rfile = uproot.open(args.infile)
if args.asimov:
    asimov = uproot.open(args.asimov)
else:
    asimov=None

input_subdir = args.infile.split("/")[-2]

cms_decor = "Preliminary"

def get_pulls(rtfile, uncertainties=None, poi_type="mu"):
    logger.info(f"get pulls and constraints")

    fitresult = rtfile["fitresults"]

    histname = f"nuisance_impact_{poi_type}"

    impacts = rtfile[histname].to_hist()

    logger.debug(f"Load pulls & constraints")
    # pick pulls & constraints
    if uncertainties is None:
        names = np.array([n for n in impacts.axes[1]])
        pulls = np.array([fitresult[n].array()[0] for n in impacts.axes[1]])
        constraints = np.array([fitresult[n+"_err"].array()[0] for n in impacts.axes[1]])
    else:
        names = np.array([n for n in impacts.axes[1]])
        pulls = np.array([fitresult[n].array()[0] for n in impacts.axes[1] if k in uncertainties])
        constraints = np.array([fitresult[n+"_err"].array()[0] for n in impacts.axes[1] if k in uncertainties])

    return names, pulls, constraints

translate = {}

def plot_pulls(rtfile, asmiov=None, max_nuisances=50):
    names, pulls, constraints = get_pulls(rfile)

    if asmiov != None:
        a_names, a_pulls, a_constraints = get_pulls(asmiov)

    other_indices = np.array([False]*len(names))
    for g, f in [
        # ("general", ["CMS_Top", "CMS_VV", "lumi", "massShift20MeV"]),
        ("scale", lambda x: x.startswith("CMS_scale")),
        ("prefire", lambda x: x.startswith("CMS_prefire")),
        ("recoil", lambda x: x.startswith("recoil")),
        ("pdf", lambda x: x.startswith("pdf")),
        ("scetlib", lambda x: x.startswith("scetlib")),
        ("resum", lambda x: x.startswith("resum")),
        ("angularCoefficients", lambda x: "AngCoeff" in x),
        ("eff_stat", lambda x: x.startswith("effStat")),
        ("eff_syst", lambda x: x.startswith("effSyst")),
        ("others", None)
    ]:
        if args.groups is not None and g not in args.groups:
            continue

        if g == "others":
            indices = ~other_indices
        elif isinstance(f, list):
            indices = np.array([n in f for n in names])
        else:
            indices = np.array([f(n) for n in names])
        
        other_indices = other_indices | indices

        g_names = names[indices]
        g_pulls = pulls[indices]
        g_constraints = constraints[indices]

        if asmiov != None:
            g_a_names = a_names[indices]
            g_a_pulls = a_pulls[indices]
            g_a_constraints = a_constraints[indices]

        n = len(g_names)
        if n <= 0:
            logger.warning(f"No match found for {g}! Continue with next one.")
            continue
        else:
            logger.debug(f"Make pull plot for {g}")

        for ni in range(int(n/max_nuisances)+1):
            first = max_nuisances * ni
            last  = min(max_nuisances * (ni+1), n)
            
            n_nuisances = last-first
            if n_nuisances <= 0:
                continue

            i_names = g_names[first: last]
            i_pulls = g_pulls[first: last]
            i_constraints = g_constraints[first: last]

            if asmiov != None:
                i_a_names = g_a_names[first: last]
                i_a_pulls = g_a_pulls[first: last]
                i_a_constraints = g_a_constraints[first: last]

            y = np.arange(n_nuisances)
            x = i_pulls
            x_err = i_constraints

            max_x = max(max(x+x_err), 1)
            min_x = min(min(x-x_err), -max_x)
            max_x = max(max_x, -min_x)

            min_y = -2.0
            max_y = n_nuisances + 1.0

            plt.close()

            fig_height = 8*(3+n_nuisances)/(3+max_nuisances) + 2
            fig = plt.figure(figsize=(6.0, fig_height))
            ax1 = fig.add_subplot(111)
            fig.subplots_adjust(hspace=0.0, left=0.4, right=0.98, top=1.0, bottom=0.1 * 10/fig_height)

            if asmiov != None:
                plt.bar(i_a_pulls, bottom=y-0.4, height=0.8, width=i_constraints*2 , color="grey", alpha=0.5)

            plt.plot([-1,-1], [min_y, max_y], color="grey", linestyle="--")
            plt.plot([0,0], [min_y, max_y], color="grey", linestyle="-")
            plt.plot([1,1], [min_y, max_y], color="grey", linestyle="--")

            plt.errorbar(x, y, xerr=x_err, color="black", linestyle='', marker=".", capsize=1.0)

            ax1.set_yticks(y)
            ax1.set_yticklabels(i_names, fontsize=12)

            range_x = max_x - min_x
            ax1.set_xlim([min_x-range_x*0.1, max_x+range_x*0.1])
            ax1.set_ylim([min_y, max_y])

            ax1.yaxis.set_minor_locator(ticker.NullLocator())

            ax1.set_xlabel("Pulls")

            outfile = f"{input_subdir}_pulls_{g}_{ni}"

            outfile += (f"_{args.postfix}" if args.postfix else "")
            plot_tools.save_pdf_and_png(outdir, outfile)

            plot_tools.save_pdf_and_png(outdir, outfile)


plot_pulls(rfile, asimov)

if output_tools.is_eosuser_path(args.outpath) and args.eoscp:
    output_tools.copy_to_eos(args.outpath, args.outfolder)