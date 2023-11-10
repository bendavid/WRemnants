
import h5py
import hist
import numpy as np

from narf import ioutils
from utilities import logging

logger = logging.child_logger(__name__)

def get_fitresult(fitresult_filename):
    h5file = h5py.File(fitresult_filename, mode='r')
    h5results = ioutils.pickle_load_h5py(h5file["results"])
    return h5results

def load_histograms(fitresult, name, axes_list, flow=False):
    logger.debug(f"Load histogram {name}")

    # axes_list contains the axes for each channel
    hflat = fitresult[name].get()
    shape_flat = hflat.view(flow=flow).shape
    nprocs = 1 if len(shape_flat)==1 else shape_flat[1]
    hnews_stack = []
    for i in range(nprocs):
        ibin=0
        hnews = []
        for axes in axes_list:
            hnew = hist.Hist(*axes, storage=hflat.storage_type())
            start = ibin
            stop = ibin+np.product(hnew.shape)
            idxs = [slice(start,stop)]
            if nprocs > 1:
                idxs.append(i)
            nbins = stop-start
            if np.product(hnew.shape) != nbins:
                raise ValueError(f"Number of bins in flat histogram ({nbins}) does not match with required shape ({hnew.shape}) and number of bins ({np.product(hnew.shape)}).")
            hnew.view(flow=flow)[...] = np.reshape(hflat.view(flow=flow)[*idxs], hnew.shape)
            hnews.append(hnew)
            ibin += np.product(hnew.shape)
        hnews_stack.append(hnews)

    if nprocs==1:
        return hnews_stack[0]

    # make channel the outer dimension
    return [x for x in zip(*hnews_stack)]
