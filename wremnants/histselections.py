from utilities import boostHistHelpers as hh
import hist
import numpy as np
from scripts.analysisTools.plotUtils.utility import safeOpenFile, safeGetObject

hist_map = {
    "eta_pt" : "nominal",
    "eta" : "nominal",
    "pt" : "nominal",
    "mll" : "nominal",
    "ptll" : "nominal",
    "ptll_mll" : "nominal",
}

def fakeHistABCD(h):
    return hh.multiplyHists(
        hh.divideHists(h[{"passIso" : True, "passMT" : False}], 
            h[{"passIso" : False, "passMT" : False}],
                cutoff=1
            ),
                #where=h[{"passIso" : False, "passMT" : True}].values(flow=True)>1),
        h[{"passIso" : False, "passMT" : True}], 
    )

def fakeHistIsoRegion(h, scale=1.):
    #return h[{"iso" : 0.3j, "mt" : hist.rebin(10)}]*scale
    return h[{"iso" : 4}]*scale

def fakeHistIsoRegionIntGen(h, scale=1.):
    print([ax.name for ax in h.axes])
    if not "qTgen" in [ax.name for ax in h.axes]:
        return h[{"iso" : 4}]*scale
    s = hist.tag.Slicer()
    print("Slicing")
    return h[{"iso" : 0, "qTgen" : s[::hist.sum]}]

def signalHistWmass(h, charge=None, passIso=True, passMT=True):
    sel = {"passIso" : passIso, "passMT" : passMT}
    if charge in [-1, 1]:
        sel.update({"charge" : -1j if charge < 0 else 1j})
    return h[sel]

# the following are utility wrapper functions for signalHistWmass with proper region selection
def histWmass_failMT_passIso(h, charge=None):
    return signalHistWmass(h, charge, True, False)

def histWmass_failMT_failIso(h, charge=None):
    return signalHistWmass(h, charge, False, False)

def histWmass_passMT_failIso(h, charge=None):
    return signalHistWmass(h, charge, False, True)

def histWmass_passMT_passIso(h, charge=None):
    return signalHistWmass(h, charge, True, True)

# TODO: Not all hists are made with these axes
def signalHistLowPileupW(h):
    if not "qTgen" in [ax.name for ax in h.axes]:
        return h[{"iso" : 0}]
    s = hist.tag.Slicer()
    return h[{"iso" : 0, "qTgen" : s[::hist.sum]}]
    
def signalHistLowPileupZ(h):
    return h

def unrolledHist(h, obs=["pt", "eta"]):
    hproj = h.project(*obs)
    bins = np.product(hproj.axes.size)
    newh = hist.Hist(hist.axis.Integer(0, bins), storage=hproj._storage_type())
    newh[...] = np.ravel(hproj)
    return newh

def applyCorrection(h, scale=1.0, offsetCorr=0.0, corrFile=None, corrHist=None, createNew=False):
    # originally intended to apply a correction differential in eta-pt
    # corrHist is a TH3 with eta-pt-charge
    # scale is just to apply an additional constant scaling (or only that if ever needed) before the correction from corrHist
    # offsetCorr is for utility, to add to the correction from the file, e.g. if the histogram is a scaling of x (e.g. 5%) one has to multiply the input histogram h by 1+x, so offset should be 1
    boost_corr = None
    if corrFile and corrHist:
        ## TDirectory.TContext() should restore the ROOT current directory to whatever it was before a new ROOT file was opened
        ## but it doesn't work at the moment, apparently the class miss the __enter__ member and the usage with "with" fails
        #with ROOT.TDirectory.TContext():
        f = safeOpenFile(corrFile, mode="READ")
        corr = safeGetObject(f, corrHist, detach=True)
        if offsetCorr:
            offsetHist = corr.Clone("offsetHist")
            ROOT.wrem.initializeRootHistogram(offsetHist, offsetCorr)
            corr.Add(offsetHist)
        f.Close()
        boost_corr = narf.root_to_hist(corr)
    # note: in fact hh.multiplyHists already creates a new histogram
    if createNew:
        hnew = hist.Hist(*h.axes, storage=hist.storage.Weight())
        hnew.values(flow=True)[...]    = scale * h.values(flow=True)
        hnew.variances(flow=True)[...] = scale * scale * h.variances(flow=True)
        if boost_corr:
            hnew = hh.multiplyHists(hnew, boost_corr)
        return hnew
    else:
        h.values(flow=True)[...]    *= scale
        h.variances(flow=True)[...] *= (scale * scale)
        if boost_corr:
            h = hh.multiplyHists(h, boost_corr)
        return h
