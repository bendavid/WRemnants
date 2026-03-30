import argparse
import pickle

# import lz4.frame

import matplotlib as mpl
import matplotlib.pyplot as plt

from wremnants.postprocessing.datagroups.datagroups import Datagroups
from wums import logging
import numpy as np
import hist
import wums.plot_tools

mpl.rcParams["figure.dpi"] = 300


parser = argparse.ArgumentParser()
parser.add_argument("-i", "--inputFile", type=str, required=True)
parser.add_argument("--debug", action="store_true", help="Print debug output")
args = parser.parse_args()

logger = logging.setup_logger("plot_muon_response", 4 if args.debug else 3)

# nominalName = "hist_qopr"

histnames = []
histnames.append("hist_qopr")
histnames.append("hist_qopr_shifted")
# histnames.append("hist_qopr_scaled_weight")
# histnames.append("hist_qopr_scaled_weight_gaussian")
# histnames.append("hNValidPixelHitsNonTrig")
# histnames.append("hist_qopr_smeared")
# histnames.append("hist_qopr_smeared_weight")
# histnames.append("hist_qopr_smeared_weight_gaussian")
histnames.append("hist_qopr_scaled_weight_moment")

datagroups = Datagroups(args.inputFile)

# if datagroups.mode != "z_dilepton":
    # raise ValueError("Expected input is the output from the dilepton histmaker")

for histname in histnames:
    datagroups.loadHistsForDatagroups(histname, syst="")


# groups = datagroups.getDatagroups()
# results = datagroups.resultsDict()
groups = datagroups.groups
# print(groups)
# print(results["Wminusmunu_2016PostVFP"].keys())

print(groups.keys())

# quit()
hists = []

for h in histnames:
    procs = datagroups.getNames()
    # print("procs", procs)
    # for proc in procs:
        # print(groups[proc].hists.keys())
    prochists = [groups[proc].hists[h] for proc in procs]
    prochists = [h for h in prochists if h is not None]
    # print(hists)
    # print(prochists)
    h = sum(prochists)

    # h = h.project("qopr")
    # print(h)

    print(np.sum(h.values(flow=True)))

    hists.append(h)
    # for proc in procs:
        # print(groups[proc].hists)
    # h = sum([groups[proc].hists[hist] for proc in procs])
    # print(h)
    # h = sum(datagroups.getDatagroupsForHist(hist))

# quit()


# hists[-1].values(flow=True)[...] *= 5e-5**3*1./np.where(hists[0].values(flow=True)!=0., hists[0].values(flow=True), 1.)
# # hists[-1].values(flow=True)[...] *= -12.*1./5e-5**2
# hists[-1].values(flow=True)[...] -= hists[0].values(flow=True)
# hists[-1].values(flow=True)[...] *= 12
# hists[-1].values(flow=True)[...] += hists[0].values(flow=True)

hists = [h.project("qopr") for h in hists]

dr = 1e-2

plt.figure()
wums.plot_tools.makePlotWithRatioToRef(hists, labels=histnames, rrange = [[1.-dr, 1.+dr]])
plt.savefig("response.pdf")


quit()

plt.figure()

sumwnom = np.sum(hists[0].values(flow=True))
# hists = [hist*(sumwnom/np.sum(hist.values(flow=True))) for hist in hists]




for h in hists:
    print(np.sum(h.values(flow=True)))

if True:
    for h, name in zip(hists, histnames):
        # h = h[...,::hist.rebin(5)]
        h = h.project("qopr")
        # h = h[12,12,0, ...]
        # h = h[::hist.rebin(10)]
        # print(h)
        h.plot(label=name)
else:
    h = hists[-1].project("qopr")
    h.plot()

# plt.xlim(0.9, 1.1)
plt.legend()

plt.savefig("response.pdf")

# hNValidPixelHitsTrig_mc = (
#     groups["Zmumu"].hists["hNValidPixelHitsTrig"]
#     + groups["Ztautau"].hists["hNValidPixelHitsTrig"]
# )
#
# hNValidPixelHitsNonTrig_mc = (
#     groups["Zmumu"].hists["hNValidPixelHitsNonTrig"]
#     + groups["Ztautau"].hists["hNValidPixelHitsNonTrig"]
# )
#
#
# print(hNValidPixelHitsTrig_mc)
#
# res = {
#     "hNValidPixelHitsTrig_data": groups["Data"].hists["hNValidPixelHitsTrig"],
#     "hNValidPixelHitsNonTrig_data": groups["Data"].hists["hNValidPixelHitsNonTrig"],
#     "hNValidPixelHitsTrig_mc": hNValidPixelHitsTrig_mc,
#     "hNValidPixelHitsNonTrig_mc": hNValidPixelHitsNonTrig_mc,
# }
#
# with lz4.frame.open("pixelcorr.pkl.lz4", "wb") as fout:
#     pickle.dump(res, fout, protocol=pickle.HIGHEST_PROTOCOL)
