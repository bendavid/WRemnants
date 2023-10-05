import ROOT
import pathlib
import hist
import narf.clingutils
import uproot
import pathlib
import hist
import pickle
import lz4.frame
from .correctionsTensor_helper import makeCorrectionsTensor
from .theory_tools import axis_helicity, moments_to_angular_coeffs
from utilities import common, logging
from utilities import boostHistHelpers as hh
import numpy as np
import h5py
import hdf5plugin
import narf

logger = logging.child_logger(__name__)

narf.clingutils.Declare('#include "syst_helicity_utils.h"')

data_dir = f"{pathlib.Path(__file__).parent}/data/"

#UL, A0...A4
axis_helicity_multidim = hist.axis.Integer(-1, 5, name="helicitySig", overflow=False, underflow=False)

#creates the helicity weight tensor
def makehelicityWeightHelper(is_w_like = False, filename=None):
    if filename is None:
        filename = f"{common.data_dir}/angularCoefficients/w_z_coeffs_scetlib_dyturboCorr.hdf5" 
    with h5py.File(filename, "r") as ff:
        out = narf.ioutils.pickle_load_h5py(ff["results"])

    corrh = out["Z"] if is_w_like else out["W"]
    if 'muRfact' in corrh.axes.name:
        corrh = corrh[{'muRfact' : 1.j,}]
    if 'muFfact' in corrh.axes.name:
        corrh = corrh[{'muFfact' : 1.j,}]
    missing = set(["y", "ptVgen", "chargeVgen", "helicity", "massVgen"]).difference(set(corrh.axes.name))
    if missing != set():
        raise ValueError (f"Axes {missing} are not present in the coeff histogram")
    
    corrh = corrh.project('massVgen','y','ptVgen','chargeVgen', 'helicity')
    if np.count_nonzero(corrh[{"helicity" : -1.j}] == 0):
        logger.warning("Zeros in sigma UL for the angular coefficients will give undefined behaviour!")
    # histogram has to be without errors to load the tensor directly
    corrh_noerrs = hist.Hist(*corrh.axes, storage=hist.storage.Double())
    corrh_noerrs.values(flow=True)[...] = corrh.values(flow=True)

    return makeCorrectionsTensor(corrh_noerrs, ROOT.wrem.WeightByHelicityHelper, tensor_rank=1)

#Muon eff vars
def make_muon_eff_stat_helpers_helicity(helper_stat, nhelicity=6):
    axes = helper_stat.tensor_axes
    nEta = axes[0].size
    nPt = axes[1].size
    nCharge = axes[2].size
    helper_stat = ROOT.wrem.muon_eff_helper_stat_helicity[nEta, nPt, nCharge, nhelicity]()
    tensor_axes = [axis_helicity_multidim, *axes]

    return helper_stat, tensor_axes

#1D tensor
# axis_all = hist.axis.Integer(0, 5, underflow = False, overflow = False, name = "reco-tracking-idip-trigger-iso")
def make_muon_eff_syst_helper_helicity(helper_syst, nhelicity=6):
    nsize=helper_syst.tensor_axes[0].size
    nvars=helper_syst.tensor_axes[1].size
    helper_syst_helicity=ROOT.wrem.tensorRank2_helper_helicity[nsize, nvars, nhelicity]()
    tensor_axes=[axis_helicity_multidim, *helper_syst.tensor_axes]
    return helper_syst_helicity, tensor_axes

#mass weights
def make_massweight_helper_helicity(mass_axis, nhelicity=6):
    tensor_axes=[axis_helicity_multidim, mass_axis]
    helper = ROOT.wrem.tensor1D_helper_helicity[mass_axis.size, nhelicity]()
    return helper, tensor_axes

#muon prefire
#this is helcity X <up/down> 
def make_muon_prefiring_helper_syst_byHelicity(nhelicity=6):
    helper_syst = ROOT.wrem.tensor1D_helper_helicity[2, nhelicity]()
    axis_tensor = [axis_helicity_multidim, common.down_up_axis]
    return helper_syst, axis_tensor

#this is helicity X <Neta,2> type
def make_muon_prefiring_helper_stat_byHelicity(helper_stat, nhelicity=6):
    nEta = helper_stat.tensor_axes[0].size
    helper_stat_helicity = ROOT.wrem.tensorupdownvar_helper_helicity[nEta, nhelicity]()
    tensor_axes = [axis_helicity_multidim, *helper_stat.tensor_axes]
    return helper_stat_helicity, tensor_axes


#for muonscale_hist
def make_dummymuonscale_helper_helicity(nweights, netabins, haxes,nhelicity=6):
    helper = ROOT.wrem.tensorRank2_helper_helicity[nweights, netabins, nhelicity]()
    tensor_axes = [axis_helicity_multidim, *haxes]
    return helper, tensor_axes

##for pdf
def make_pdfweight_helper_helicity(npdf, pdf_axes, nhelicity=6):
    helper=ROOT.wrem.tensor1D_helper_helicity[npdf, nhelicity]()
    tensor_axes=[axis_helicity_multidim,pdf_axes]
    return helper, tensor_axes

#for qcd scale
def make_qcdscale_helper_helicity(qcd_axes,nhelicity=6):
    helper = ROOT.wrem.tensorRank2_helper_helicity[3, 3, nhelicity]()
    tensor_axes = [axis_helicity_multidim, *qcd_axes]
    return helper, tensor_axes
