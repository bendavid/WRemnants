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

if "A-M-separated" in args.input:
    A_constraints = [getattr(tree, label+'_err') for label in labels if 'non_closure_parametrized_A' in label]
    A_avg = sum(A_constraints) / len(A_constraints)
    print("the average of constraints on A uncs is: ", A_avg)

    M_constraints = [getattr(tree, label+'_err') for label in labels if 'non_closure_parametrized_M' in label]
    M_avg = sum(M_constraints) / len(M_constraints)
    print("the average of constraints on M uncs is: ", M_avg)

elif "A-M-combined" in args.input:
    param_constraints = [getattr(tree, label+'_err') for label in labels if 'non_closure_parametrized' in label]
    param_avg = sum(param_constraints) / len(param_constraints)
    print("the average of constraints on param uncs is: ", param_avg)
elif "binned" in args.input:
    binned_constraints = [getattr(tree, label+'_err') for label in labels if 'non_closure_binned' in label]
    binned_avg = sum(binned_constraints) / len(binned_constraints)
    print("the average of constraints on binned uncs is: ", binned_avg)
