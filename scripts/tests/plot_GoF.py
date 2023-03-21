# plot goodness of fit from toys with data and/or pseudodata 
# use full nll as test statistic

import uproot
import os
import argparse
import pdb
import matplotlib.pyplot as plt
import numpy as np
from scipy.stats import chi2

parser = argparse.ArgumentParser()
parser.add_argument("-i","--input", type=str, help="fitresult root file with toys")
parser.add_argument("-t","--testData", default=[], nargs="*" , help="fitresult root file with test data (data, pseudodata)")
parser.add_argument("-l","--testDabels", default=[], nargs="*" , help="labels for test data (data, pseudodata) for plotting")
parser.add_argument("-a","--asimov", type=str, default=None , help="fitresult root file with asimov data for normalization")
parser.add_argument("-n","--nDoF", default=20,type=int, help="Number of degrees of freedom, for the chi2 function in the plot")
parser.add_argument("-o","--outputDir", default="./",type=str, help="output directory for plots")
args = parser.parse_args()

if not os.path.isdir(args.outputDir):
    os.mkdir(args.outputDir)

with uproot.open(f"{args.input}:fitresults") as utree:
    nll = utree['nllvalfull'].array(library="np")
    nll_saturated = utree['satnllvalfull'].array(library="np")
    # deviance
    toys = 2*(nll - nll_saturated)

data=[]
for tdata in args.testData:
    with uproot.open(f"{tdata}:fitresults") as utree:
        nll = utree['nllvalfull'].array(library="np")[0]
        nll_saturated = utree['satnllvalfull'].array(library="np")[0]
        data.append(2*(nll - nll_saturated))
    
data = np.array(data)

if args.testLabels ==[]:
    testlabels = ["" for _ in data]
else:
    testlabels = args.testLabels

# calculate p value as fraction of toys that have a higher likelihood than data
pvalues = np.array([sum(toys>dt)/len(toys) for dt in data])

if args.asimov:
    # we subtract the asimov data to subtract the contribution from the bin-by-bin statistial MC uncertainty,
    #   for simplicity this is not included in the saturated model but is the same for all toys/ data/ pseudodata for a given model
    with uproot.open(f"{args.asimov}:fitresults") as utree:
        nll = utree['nllvalfull'].array(library="np")[0]
        nll_saturated = utree['satnllvalfull'].array(library="np")[0]
        asimov = 2*(nll - nll_saturated)
        toys -= asimov
        data -= asimov


xlabel="$2 (\\mathrm{log}(\\mathcal{L}_\\mathrm{S}) - \\mathrm{log}(\\mathcal{L}_\\mathrm{P})) $"
colors = ["red", "blue", "green"]
saveas="png"
textsize = 16
labelsize = 12
nBins = 50

rangex = min(np.append(toys,data)), max(np.append(toys,data))
rangex = rangex[0]-(0.01*(rangex[1]-rangex[0])), rangex[1]+(0.01*(rangex[1]-rangex[0])) 

plt.clf()
fig = plt.figure()
fig.subplots_adjust(left=0.15, right=0.99, top=0.99, bottom=0.125)
ax = fig.add_subplot(111)

nEntries, bins, _ = ax.hist(toys, bins=nBins, range=rangex, color="yellow")
ax.step(bins, np.append(nEntries,nEntries[-1]), color="black", where="post", linewidth=0.5)

rangey = 0, max(nEntries)*1.2

for i, dt, l, p in zip(range(len(data)), data, testlabels, pvalues):
    ax.plot([dt,dt], rangey, color=colors[i], label=f"{l} $p={p*100}\%$")

norm = (rangex[1]-rangex[0])/nBins*len(toys)

x = np.linspace(*rangex, 100)
# ax.plot(x, chi2.pdf(x, args.nDoF)*norm, 'k-', lw=1, label='chi2 pdf')

ax.set_ylabel("Number of toys", fontsize=textsize)

ax.set_xlabel(xlabel, fontsize=textsize)

ax.set_xlim(rangex)
ax.set_ylim(rangey)

plt.legend()

plt.xticks(fontsize = labelsize)
plt.yticks(fontsize = labelsize)

histname = f"{args.outputDir}/hist_GoF.{saveas}"
print(f"Save: {histname}")
plt.savefig(histname)
plt.close()
