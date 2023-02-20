import setupCombineWMass,setupLowPU_Z

parser = setupCombineWMass.make_parser()
parser = setupLowPU_Z.make_parser(parser)
parser.add_argument("--inputFileLowPU", type=str, default="", help="Input file for low pu step")
parser = common.set_parser_default(parser, "baseDir", "combineResults/WMassLowPUComb")

args = parser.parse_args()

base_out = args.outfolder
lowpu_out = "/".join([base_out, "lowPU"])
highpu_out = "/".join([base_out, "highPU"])
args.outfolder = highpu_out
setupCombineWMass.main(args)

args.inputFile = args.inputFileLowPU
args.outfolder = lowpu_out
# TODO: Call W if appropriate
if args.wlike:
    setupLowPU_Z.main(args)

with open("/".join([args.baseDir, base_out, "combine.sh"]), "w") as f:
    f.write("combineCards.py --noDirPrefix ")
    f.write(f"lowPU=lowPU/lowPU_Z{args.flavor}_{args.met}_wmass.txt ")
    f.write(f"minus=highPU/ZMassWLike_minus.txt ")
    f.write(f"plus=highPU/ZMassWLike_plus.txt > ZMassWLike_lowPU_comb.txt\n")
    f.write("text2hdf5.py --X-allow-no-signal ZMassWLike_lowPU_comb.txt\n")
    f.write("combinetf.py -t -1 --doImpacts --binByBinStat ZMassWLike_lowPU_comb.hdf5")
