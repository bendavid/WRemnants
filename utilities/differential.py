from utilities import logging, common
import hist

logger = logging.child_logger(__name__)

axis_xnorm = hist.axis.Integer(0, 1, underflow=False, overflow=False, name = "count")
axis_fiducial = hist.axis.Boolean(name = "fiducial")

eta_binning = [0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.3, 1.5, 1.7, 1.9, 2.1, 2.4] # 18 eta bins

def get_pt_eta_axes(n_bins_pt, min_pt, max_pt, n_bins_eta=0):

    # gen axes for differential measurement
    if n_bins_eta > 0:
        axis_etaGen = hist.axis.Regular(n_bins_eta, 0, 2.4, underflow=False, overflow=False, name = "etaGen")
    else:
        axis_etaGen = hist.axis.Variable(eta_binning, underflow=False, overflow=False, name = "etaGen")
    axis_ptGen = hist.axis.Regular(n_bins_pt, min_pt, max_pt, underflow=False, overflow=False, name = "ptGen")    

    axes = [axis_ptGen, axis_etaGen, axis_fiducial]
    cols = ["ptGen", "etaGen", "fiducial"]

    return axes, cols

def get_pt_eta_charge_axes(n_bins_pt, min_pt, max_pt, n_bins_eta=0):

    axes, cols = get_pt_eta_axes(n_bins_pt, min_pt, max_pt, n_bins_eta)

    axis_qGen = hist.axis.Regular(2, -2., 2., underflow=False, overflow=False, name = "qGen")

    axes.append(axis_qGen)
    cols.append("qGen")

    return axes, cols

def get_ptV_axes(binning):
    # axes for fiducial measurement of Z in dilepton pT(l,l)
    axis_ptVGen = hist.axis.Variable(common.ptV_binning, underflow=False, overflow=False, name = "ptVGen")

    axes = [axis_ptVGen, axis_fiducial]
    cols = ["ptVGen", "fiducial"]

    return axes, cols