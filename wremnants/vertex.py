import ROOT
import pathlib
import hist
import narf
import numpy as np
import boost_histogram as bh
from utilities import common
from utilities import common, logging

narf.clingutils.Declare('#include "vertex.h"')

logger = logging.child_logger(__name__)

data_dir = common.data_dir

def make_vertex_helper(era = None, filename = None, dataYear = 2016):

    eradict = { "2016PreVFP" :  "BtoF",
                "2016PostVFP" : "GtoH",
    }
    filedict = {
        2016 : "/vertex/vertexPileupWeights.root",
        2017 : "/vertex/vertexPileupWeights_2017.root",
        2018 : "/vertex/vertexPileupWeights_2018.root"
    }

    hnamedict = {
        2016 : f"weight_vertexZ_pileup_{eradict[era]}",
        2017 : f"weight_vertexZ_pileup_{dataYear}",
        2018 : f"weight_vertexZ_pileup_{dataYear}"
    }

    if filename is None:
        filename = data_dir + filedict[dataYear]
        print("Vertex weight fname:", filename)
    logger.debug(f"vertex.py: will read weight_vertexZ_pileup_{eradict[era]} from {filename}")
    fmc = ROOT.TFile.Open(filename)
    mchist = fmc.Get(hnamedict[dataYear])
    mchist.SetDirectory(0)
    fmc.Close()

    ## for the histogram of preVFP, last PU bin in [40-45] is empty because of very small stat
    ## better to fill that bin with content of previous one, otherwise we are effectively cutting events based on PU
    if era == "2016PreVFP":
        lastPUbin = mchist.GetNbinsY()
        for ix in range(1, mchist.GetNbinsX() +1):
            mchist.SetBinContent(ix, lastPUbin, mchist.GetBinContent(ix, lastPUbin-1))
            mchist.SetBinError(  ix, lastPUbin, mchist.GetBinError(  ix, lastPUbin-1))
            
    helper = ROOT.wrem.vertex_helper(mchist)

    return helper
