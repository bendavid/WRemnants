from utilities import logging
import hist

logger = logging.child_logger(__name__)

axis_xnorm = hist.axis.Regular(1, 0., 1., name = "count", underflow=False, overflow=False)

def get_pt_eta_axes(n_bins, min_pt, max_pt, max_eta):

    # gen axes for differential measurement
    if len(n_bins) == 2:
        genBinsPt = n_bins[0]
        genBinsEta = n_bins[1]
    else:
        raise IOError(f"Specified format '--genBins {args.genBins}' can not be processed! Please specify the number of gen bins for pT and |eta|")
    axis_ptGen = hist.axis.Regular(genBinsPt, min_pt, max_pt, name = "ptGen")
    axis_etaGen = hist.axis.Regular(genBinsEta, 0, max_eta, underflow=False, name = "etaGen")
    axes = [axis_ptGen, axis_etaGen]
    cols = ["ptGen", "etaGen"]

    return axes, cols

def get_pt_eta_charge_axes(n_bins, min_pt, max_pt, max_eta):

    axes, cols = get_pt_eta_axes(n_bins, min_pt, max_pt, max_eta)

    axis_qGen = hist.axis.Regular(2, -2., 2., underflow=False, overflow=False, name = "qGen")

    axes.append(axis_qGen)
    cols.append("qGen")

    return axes, cols
