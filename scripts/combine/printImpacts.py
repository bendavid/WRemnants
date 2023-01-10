import uproot
import argparse
from utilities import input_tools

def parseArgs():
    parser = argparse.ArgumentParser()
    parser.add_argument("-u", "--ungroup", action='store_true', help="Use ungrouped nuisances")
    parser.add_argument("-n", "--nuisance", type=str, help="Only print value for specific nuiance")
    parser.add_argument("-s", "--sort", action='store_true', help="Sort nuisances by impact")
    parser.add_argument("-f", "--inputFile", 
        default="/Users/kenneth/cernbox/CombineStudies/WGen/etal_ptl_smear_unrolled_scetlib/fitresults_123456789.root", 
        help="fitresults output ROOT file from combinetf")
    return parser.parse_args()

args = parseArgs()
rtfile = uproot.open(args.inputFile)
impacts,labels = input_tools.readImpacts(rtfile, not args.ungroup, sort=args.sort)
if args.nuisance:
    if args.nuisance not in labels:
        raise ValueError(f"Invalid nuisance {args.nuisance}. Options are {labels}")
    print(f"Impact of nuisance {args.nuisance} is {impacts[labels.index(args.nuisance)]} MeV")
else:
    print("Impact of all systematics (in MeV)")
    print("\n".join([f"   {k}: {round(v*100, 2)}" for k,v in zip(labels, impacts)]))
