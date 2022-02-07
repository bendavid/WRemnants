import ROOT
import pathlib
import hist
import narf
import numpy as np

ROOT.gInterpreter.Declare('#include "theory_corrections.h"')

data_dir = f"{pathlib.Path(__file__).parent}/data/"

def makeScetlibCorrHelper(filename=f"{data_dir}/N3LLCorrections/inclusive_Wp_pT.npz"):
    corrf = np.load(filename, allow_pickle=True)
    bins = corrf["bins"]
    axis_syst = hist.axis.Regular(len(bins[0]) - 1, bins[0][0], bins[0][-1], 
                    name="systIdx", overflow=False, underflow=False)
    axis_mass = hist.axis.Variable(bins[1], name="mass", overflow=False, underflow=False)
    axis_y = hist.axis.Variable(bins[2], name="y", overflow=False, underflow=False)
    axis_pt = hist.axis.Regular(len(bins[-1]) - 1, bins[-1][0], bins[-1][-1], name="pT", 
                overflow=False, underflow=False)

    corrh = hist.Hist(axis_mass, axis_y, axis_pt, axis_syst)
    corrh[...] = np.moveaxis(corrf["scetlibCorr3D_Wp"], 0, -1)

    corrhConv = narf.hist_to_pyroot_boost(corrh, tensor_rank=1)
    helper = ROOT.wrem.ScetlibCorrectionsHelper[type(corrhConv).__cpp_name__](ROOT.std.move(corrhConv))
    return helper
