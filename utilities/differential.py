from utilities import logging, common
import hist

logger = logging.child_logger(__name__)

eta_binning = [0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.3, 1.5, 1.7, 1.9, 2.1, 2.4] # 18 eta bins

def get_pt_eta_axes(n_bins_pt, min_pt, max_pt, n_bins_eta=0):

    # gen axes for differential measurement
    if n_bins_eta > 0:
        axis_absEtaGen = hist.axis.Regular(n_bins_eta, 0, 2.4, underflow=False, overflow=False, name = "absEtaGen")
    else:
        axis_absEtaGen = hist.axis.Variable(eta_binning, underflow=False, overflow=False, name = "absEtaGen")
    axis_ptGen = hist.axis.Regular(n_bins_pt, min_pt, max_pt, underflow=False, overflow=False, name = "ptGen")    

    axes = [axis_ptGen, axis_absEtaGen]
    cols = ["ptGen", "absEtaGen"]

    return axes, cols

def get_pt_eta_charge_axes(n_bins_pt, min_pt, max_pt, n_bins_eta=0):

    axes, cols = get_pt_eta_axes(n_bins_pt, min_pt, max_pt, n_bins_eta)

    axis_qGen = hist.axis.Regular(2, -2., 2., underflow=False, overflow=False, name = "qGen")

    axes.append(axis_qGen)
    cols.append("qGen")

    return axes, cols

def get_dilepton_axes(gen_vars, gen_axes):
    # axes for fiducial measurement of Z in dilepton e.g. pT(Z), |yZ|

    axes = []
    cols = []
    selections = []

    for var in gen_vars:
        axes.append(gen_axes[var])
        cols.append(var)
    
    # selections for out of fiducial region
    if "ptVGen" in gen_vars:
        selections.append("ptVGen < {0}".format(gen_axes["ptVGen"].edges[-1]))

    if "absYVGen" in gen_vars:
        selections.append("absYVGen < {0}".format(gen_axes["absYVGen"].edges[-1]))

    return axes, cols, selections