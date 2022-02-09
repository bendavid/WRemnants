import ROOT
import pathlib

ROOT.gInterpreter.AddIncludePath(f"{pathlib.Path(__file__).parent}/include/")

ROOT.gInterpreter.Declare('#include "muonCorr.h"')
ROOT.gInterpreter.Declare('#include "theoryTools.h"')
ROOT.gInterpreter.Declare('#include "pileupWeights.h"')
ROOT.gInterpreter.Declare('#include "utils.h"')

from .datasets import datasets2016

from .muon_prefiring import make_muon_prefiring_helpers
from .muon_efficiencies import make_muon_efficiency_helpers
from .scetlib_corrections import makeScetlibCorrHelper

data_dir = f"{pathlib.Path(__file__).parent}/data/"
