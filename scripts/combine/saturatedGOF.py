import ROOT
import scipy
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("infile", type=str, help="combinetf output ROOT file")
args = parser.parse_args()

rtfile = ROOT.TFile(args.infile)
tree = rtfile.Get("fitresults")
tree.GetEntry(0)
nbins = rtfile.Get("obs").GetNbinsX()

print(f"-loglikelihood_{{full}} = {tree.nllvalfull}")
print(f"-loglikelihood_{{saturated}} = {tree.satnllvalfull}")
print(f"2*(nllfull-nllsat) = {2*(tree.nllvalfull-tree.satnllvalfull)}")
print(f"nbins = {nbins}")
print("chi2 probability =", scipy.stats.chi2.sf(2*(tree.nllvalfull-tree.satnllvalfull), nbins))
