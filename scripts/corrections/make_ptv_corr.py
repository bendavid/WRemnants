from wremnants.datasets.datagroups import Datagroups
from utilities import boostHistHelpers as hh, common, logging
from utilities.io_tools import input_tools, output_tools
import numpy as np
import hist
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("-i", "--inputFile", type=str, required=True)
parser.add_argument("--outpath", type=str, default=f"{common.data_dir}/TheoryCorrections", help="Output path")
parser.add_argument("--debug", action='store_true', help="Print debug output")
parser.add_argument("-p", "--postfix", type=str, help="Postfix for output file name", default=None)
parser.add_argument("--proc", type=str, required=True, choices=["z", "w", ], help="Process")
args = parser.parse_args()

logger = logging.setup_logger("make_ptv_corr", 4 if args.debug else 3)

nominalName = "nominal"
syst = "uncorr"
hist_name = syst
datagroups = Datagroups(args.inputFile)

if datagroups.mode != "dilepton":
    raise ValueError("Expected input is the output from the dilepton histmaker")

datagroups.setNominalName(nominalName)
datagroups.loadHistsForDatagroups(nominalName, syst=syst)
datagroups.addSummedProc(syst, name=syst, rename="MC_sum")

groups = datagroups.getDatagroups()
ratio = hh.divideHists(*[groups[x].hists[hist_name].project("ptll") for x in ['Data', 'MC_sum']])

pt_ax = ratio.axes["ptll"]
pt_ax._ax.metadata['name'] = "ptVgen"
charge_ax = hist.axis.Regular(*(1, -1, 1) if args.proc == 'z' else (2, -2, 2), name="chargeVgen", flow=False)

corrh = hist.Hist(hist.axis.Regular(1, 0, 13000, name="massVgen", flow=False), 
                  hist.axis.Regular(1, 0, 10, name="absYVgen", underflow=False, overflow=True), 
                  pt_ax, charge_ax, hist.axis.Regular(1, 0, 1, name="vars"))

corrh[...] = ratio.values(flow=True)[np.newaxis,np.newaxis,:,np.newaxis, np.newaxis]

output_dict = {
    "MC_data_ratio" : corrh,
    "data_hist" : groups['Data'].hists[hist_name],
    "MC_sum_hist" : groups['MC_sum'].hists[hist_name],
}

meta_dict = input_tools.get_metadata(args.inputFile)
outfile = args.outpath+"/dataPtll"
output_tools.write_theory_corr_hist(outfile, args.proc.upper(), output_dict, args, meta_dict)
logger.info(f"Average correction is {np.mean(corrh.values())}")
