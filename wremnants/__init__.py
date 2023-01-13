import ROOT
import pathlib

ROOT.gInterpreter.AddIncludePath(f"{pathlib.Path(__file__).parent}/include/")

ROOT.gInterpreter.Declare('#include "muonCorr.h"')
ROOT.gInterpreter.Declare('#include "histoScaling.h"')
ROOT.gInterpreter.Declare('#include "histHelpers.h"')
ROOT.gInterpreter.Declare('#include "utils.h"')
ROOT.gInterpreter.Declare('#include "csVariables.h"')
ROOT.gInterpreter.Declare('#include "EtaPtCorrelatedEfficiency.h"')

from .datasets import datasets2016
from .datasets import datasetsLowPU

from .muon_prefiring import make_muon_prefiring_helpers
from .muon_efficiencies_smooth import make_muon_efficiency_helpers_smooth
from .muon_efficiencies_binned import make_muon_efficiency_helpers_binned
from .qcdScaleByHelicity_helper import makeQCDScaleByHelicityHelper
from .pileup import make_pileup_helper
from .vertex import make_vertex_helper
from .syst_tools import scale_helicity_hist_to_variations
from .theory_tools import axis_helicity, scale_tensor_axes, define_prefsr_vars, define_ew_vars, moments_to_angular_coeffs
from .muon_calibration import *
from .muon_validation import *

data_dir = f"{pathlib.Path(__file__).parent}/data/"
