from utilities import boostHistHelpers as hh
import hist
import numpy as np

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

def signalHistWmass(h, charge=None, passIso=1, passMT=1):
    sel = {"passIso" : passIso, "passMT" : passMT}
    if charge in [-1, 1]:
        sel.update({"charge" : -1j if charge < 0 else 1j})
    return h[sel]

# the following are utility wrapper functions for signalHistWmass with proper region selection
def histWmass_failMT_passIso(h, charge=None):
    return signalHistWmass(h, charge, 1, 0)

def histWmass_failMT_failIso(h, charge=None):
    return signalHistWmass(h, charge, 0, 0)

def histWmass_passMT_failIso(h, charge=None):
    return signalHistWmass(h, charge, 0, 1)

def histWmass_passMT_passIso(h, charge=None):
    return signalHistWmass(h, charge, 1, 1)

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
    bins = np.multiply(*hproj.axes.size)
    newh = hist.Hist(hist.axis.Integer(0, bins), storage=hproj._storage_type())
    newh[...] = np.ravel(hproj)
    return newh
