import ROOT
import pathlib

ROOT.gInterpreter.AddIncludePath(f"{pathlib.Path(__file__).parent}/include/")

ROOT.gInterpreter.Declare('#include "functions.h"')
ROOT.gInterpreter.Declare('#include "functionsWMass.h"')
ROOT.gInterpreter.Declare('#include "pileupWeights.h"')

from .datasets import datasets2016

data_dir = f"{pathlib.Path(__file__).parent}/data/"
