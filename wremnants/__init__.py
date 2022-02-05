import ROOT
import pathlib

ROOT.gInterpreter.AddIncludePath(f"{pathlib.Path(__file__).parent}/include/")

ROOT.gInterpreter.Declare('#include "functions.h"')
ROOT.gInterpreter.Declare('#include "functionsWMass.h"')
ROOT.gInterpreter.Declare('#include "pileupWeights.h"')
ROOT.gInterpreter.Declare('#include "utils.h"')

from .datasets import datasets2016

from .muon_prefiring import make_muon_prefiring_helpers

data_dir = f"{pathlib.Path(__file__).parent}/data/"
