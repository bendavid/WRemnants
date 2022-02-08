import hist
import boost_histogram as bh
import numpy as np

def valsAndVariances(h1, h2, allowBroadcast):
    if not allowBroadcast or len(h1.axes) == len(h2.axes):
        return (h1.values(), h2.values(), h1.variances(), h2.variances())
    # Should probably check that all other axes are the same
    elif allowBroadcast and len(h1.axes) == len(h2.axes)-1 and h2.axes[-1].name == "systIdx":
        return (h1.values()[...,np.newaxis], h2.values(), h1.variances()[...,np.newaxis], h2.variances())
    elif allowBroadcast and len(h1.axes)-1 == len(h2.axes) and h1.axes[-1].name == "systIdx":
        return (h1.values(), h2.values()[...,np.newaxis], h1.variances(), h2.variances()[...,np.newaxis])
    else:
        raise ValueError("Incompatible hists for math operation")

def broadcastOutHist(h1, h2):
    if len(h1.axes) == len(h2.axes):
        return h1 if h1.axes[-1].size > h2.axes[-1].size else h2
    return h1 if len(h1.axes) > len(h2.axes) else h2

def divideHists(h1, h2, where=True, allowBroadcast=True):
    h1vals,h2vals,h1vars,h2vars = valsAndVariances(h1, h2, allowBroadcast)
   
    # To get the broadcast shape right
    outh = h1 if not allowBroadcast else broadcastOutHist(h1, h2)
    val = np.divide(h1vals, h2vals, out=np.zeros_like(outh), where=where)
    relvars = relVariances(h1vals, h2vals, h1vars, h2vars)
    var = val*sum(relVariances(h1vals, h2vals, h1vars, h2vars))
    var *= val
    newh = hist.Hist(*outh.axes, storage=hist.storage.Weight())
    newh[...] = np.stack((val, var), axis=-1)
    return newh

def relVariances(h1vals, h2vals, h1vars, h2vars):
    h1nonzero = np.abs(h1vals) > 1e-3
    h2nonzero = np.abs(h2vals) > 1e-3
    rel1 = np.divide(np.divide(h1vars, h1vals, out=np.zeros_like(h1vals), where=h1nonzero),
            h1vals, out=np.zeros_like(h1vals), where=h1nonzero)
    rel2 = np.divide(np.divide(h2vars, h2vals, out=np.zeros_like(h2vals), where=h2nonzero), 
            h2vals, out=np.zeros_like(h2vals), where=h2nonzero)
    return (rel1, rel2)

def multiplyHists(h1, h2, allowBroadcast=True):
    h1vals,h2vals,h1vars,h2vars = valsAndVariances(h1, h2, allowBroadcast)
    val = h1vals*h2vals
    var = val*val*sum(relVariances(h1vals, h2vals, h1vars, h2vars))

    outh = h1 if not allowBroadcast else broadcastOutHist(h1, h2)
    newh = hist.Hist(*outh.axes, storage=hist.storage.Weight())
    newh[...] = np.stack((val, var), axis=-1)
    return newh

def addHists(h1, h2, allowBroadcast=True):
    h1vals,h2vals,h1vars,h2vars = valsAndVariances(h1, h2, allowBroadcast)
    outh = h1 if not allowBroadcast else broadcastOutHist(h1, h2)

    newh = hist.Hist(*outh.axes, storage=hist.storage.Weight())
    newh[...] = np.stack((h1vals+h2vals, h1vars+h2vars), axis=-1)
    return newh

def mirrorHist(hvar, hnom, cutoff=0.1):
    hnew = multiplyHists(divideHists(hnom, hvar, where=np.abs(hvar.values())>cutoff), hnom)
    return hnew

def extendHistByMirror(hvar, hnom):
    axes = hvar.axes
    if axes[-1].name != "systIdx":
        hvar = addSystAxis(hvar)
    mirror = mirrorHist(hvar, hnom)
    offset = axes[-1].edges[0]
    ax = hist.axis.Regular(2*axes[-1].size, offset, 2*axes[-1].size+offset, name=axes[-1].name)
    hnew = hist.Hist(*axes[:-1], ax, storage=hist.storage.Weight())
    vals = np.concatenate((hvar.values(), mirror.values()), axis=-1)
    varis = np.concatenate((hvar.variances(), mirror.variances()), axis=-1) 
    hnew[...] = np.stack((vals, varis), axis=-1)
    return hnew

def addSystAxis(h, size=1, offset=0):
    hnew = hist.Hist(*h.axes,hist.axis.Regular(size,offset,size+offset, name="systIdx"), storage=hist.storage.Weight())
    # Broadcast to new shape
    newvals = hnew.values()+h.values()[...,np.newaxis]
    newvars = hnew.variances()+h.variances()[...,np.newaxis]
    hnew[...] = np.stack((newvals, newvars), axis=-1)
    return hnew

