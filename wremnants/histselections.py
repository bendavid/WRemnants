from . import boostHistHelpers as hh
import hist
import numpy as np

def fakeHistABCD(h):
    return hh.multiplyHists(
        hh.divideHists(h[{"passIso" : True, "passMT" : False}], 
            h[{"passIso" : False, "passMT" : True}],
                cutoff=1
            ),
                #where=h[{"passIso" : False, "passMT" : True}].values(flow=True)>1),
        h[{"passIso" : False, "passMT" : False}], 
    )

def fakeHistIsoRegion(h, scale=1.):
    return h[{"iso" : 0.3j, "mt" : hist.rebin(10)}]*scale

def signalHistWmass(h, charge=None):
    sel = {"passIso" : 1, "passMT" : 1}
    if charge in [-1, 1]:
        sel.update({"charge" : -1j if charge < 0 else 1j})
    return h[sel]

# TODO: Not all hists are made with these axes
def signalHistLowPileupW(h):
    return h[{"iso" : 0.j, "mt" : hist.rebin(10)}]

def unrolledHist(h, obs=["pt", "eta"]):
    bins = np.multiply(*[a.size for a in h.axes[:2]])
    newh = hist.Hist(hist.axis.Regular(bins, 0, bins), storage=hist.storage.Weight())
    newh[...] = np.ravel(h.project(*obs))
    return newh
