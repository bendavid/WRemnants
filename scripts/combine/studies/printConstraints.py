from utilities import input_tools
import ROOT
import uproot
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("-i", "--input", type=str, help="path to input fit results root file")
args = parser.parse_args()

urf = uproot.open(args.input)
rf = ROOT.TFile.Open(args.input)
impacts, labels, norm = input_tools.readImpacts(urf, None, normalize=False)
tree = rf.Get('fitresults')
tree.GetEntry(0)

constraint_dic = {label:getattr(tree, label+'_err') for label in labels if 'non_closure' in label}
for kv in (sorted(constraint_dic.items(), key = lambda kv: kv[1])): print(kv)

