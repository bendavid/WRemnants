import hist
import pathlib
import argparse
import logging

wremnants_dir = f"{pathlib.Path(__file__).parent}/../wremnants"
data_dir = f"{wremnants_dir}/data/"

wprocs = ["WplusmunuPostVFP", "WminusmunuPostVFP", "WminustaunuPostVFP", "WplustaunuPostVFP"]
zprocs = ["ZmumuPostVFP", "ZtautauPostVFP"]
vprocs = wprocs+zprocs

wprocs_lowpu = ["WminusJetsToMuNu", "WminusJetsToENu", "WminusJetsToTauNu", "WplusJetsToMuNu", "WplusJetsToENu", "WplusJetsToTauNu"]
zprocs_lowpu = ["Zmumu", "Zee", "Ztautau"]
zprocs_recoil = ["Zmumu", "Zee"]
vprocs_lowpu = wprocs_lowpu+zprocs_lowpu

# unfolding axes for low pu
axis_recoil_reco_ptZ = hist.axis.Variable([0, 5, 10, 15, 20, 30, 40, 50, 60, 75, 90, 150], name = "recoil_reco", underflow=False, overflow=True)
axis_recoil_gen_ptZ = hist.axis.Variable([0.0, 10.0, 20.0, 40.0, 60.0, 90.0, 150], name = "recoil_gen", underflow=False, overflow=True)

# standard regular axes
axis_eta = hist.axis.Regular(48, -2.4, 2.4, name = "eta")
axis_pt = hist.axis.Regular(29, 26., 55., name = "pt")
#ptV_binning = [0, 2, 3, 4, 4.75, 5.5, 6.5, 8, 9, 10, 12, 14, 16, 18, 20, 23, 27, 32, 40, 55, 100]
ptV_binning = [0, 2, 3, 4, 5, 6, 7, 8, 9, 10, 12, 14, 16, 18, 20, 23, 27, 32, 40, 55, 100]
## 5% quantiles from aMC@NLO used in SMP-18-012
#ptV_10quantiles_binning = [0.0, 1.971, 2.949, 3.838, 4.733, 5.674, 6.684, 7.781, 8.979, 10.303, 11.777, 13.435, 15.332, 17.525, 20.115, 23.245, 27.173, 32.414, 40.151, 53.858, 13000.0]
## 10% quantiles from aMC@NLO used in SMP-18-012 with some rounding <== This one worked fine with toys
ptV_10quantiles_binning = [0.0, 2.95, 4.73, 6.68, 8.98, 11.78, 15.33, 20.11, 27.17, 40.15, 120]
absYV_binning = [0, 0.25, 0.5, 0.75, 1, 1.25, 1.5, 1.75, 2, 2.25, 2.5, 2.75, 3, 3.25, 3.5, 3.75, 4]

# categorical axes in python bindings always have an overflow bin, so use a regular
# axis for the charge
axis_charge = hist.axis.Regular(2, -2., 2., underflow=False, overflow=False, name = "charge")

down_up_axis = hist.axis.Regular(2, -2., 2., underflow=False, overflow=False, name = "downUpVar")
down_nom_up_axis = hist.axis.Regular(3, -1.5, 1.5, underflow=False, overflow=False, name = "downNomUpVar")

axis_passIso = hist.axis.Boolean(name = "passIso")
axis_passMT = hist.axis.Boolean(name = "passMT")

nominal_axes = [axis_eta, axis_pt, axis_charge, axis_passIso, axis_passMT]

def common_parser():
    parser = argparse.ArgumentParser()
    parser.add_argument("-j", "--nThreads", type=int, help="number of threads")
    parser.add_argument("--debug", action='store_true', help="Debug output")
    initargs,_ = parser.parse_known_args()

    logging.basicConfig(level=logging.INFO if not initargs.debug else logging.DEBUG)

    import ROOT
    ROOT.gInterpreter.ProcessLine(".O3")
    if not initargs.nThreads:
        ROOT.ROOT.EnableImplicitMT()
    elif initargs.nThreads != 1:
        ROOT.ROOT.EnableImplicitMT(initargs.nThreads)
    import narf
    import wremnants
    from wremnants import theory_tools

    parser.add_argument("--pdfs", type=str, nargs="*", default=["nnpdf31"], choices=theory_tools.pdfMapExtended.keys(), help="PDF sets to produce error hists for")
    parser.add_argument("--altPdfOnlyCentral", action='store_true', help="Only store central value for alternate PDF sets")
    parser.add_argument("--maxFiles", type=int, help="Max number of files (per dataset)", default=-1)
    parser.add_argument("--filterProcs", type=str, nargs="*", help="Only run over processes matched by (subset) of name", default=[])
    parser.add_argument("--v8", action='store_true', help="Use NanoAODv8. Default is v9")
    parser.add_argument("-p", "--postfix", type=str, help="Postfix for output file name", default=None)
    parser.add_argument("--theory_corr", nargs="*", choices=["scetlib", "scetlibVars", "scetlibMSHT20", "scetlibHelicity", "dyturbo", "matrix_radish"], 
        help="Apply corrections from indicated generator. First will be nominal correction.", default=[])
    parser.add_argument("--theory_corr_alt_only", action='store_true', help="Save hist for correction hists but don't modify central weight")
    parser.add_argument("--skipHelicity", action='store_true', help="Skip the qcdScaleByHelicity histogram (it can be huge)")
    parser.add_argument("--eta", nargs=3, type=float, help="Eta binning as 'nbins min max' (only uniform for now)", default=[48,-2.4,2.4])
    parser.add_argument("--pt", nargs=3, type=float, help="Pt binning as 'nbins,min,max' (only uniform for now)", default=[29,26.,55.])
    parser.add_argument("--no_recoil", action='store_true', help="Don't apply recoild correction")
    parser.add_argument("--no-vertex_weight", dest="vertex_weight", action='store_false', help="Do not apply reweighting of vertex z distribution in MC to match data")
    parser.add_argument("--validationHists", action='store_true', help="make histograms used only for validations")
    return parser,initargs

def common_parser_combine():
    from wremnants import theory_tools
    parser = argparse.ArgumentParser()
    parser.add_argument("-d", "--baseDir", type=str, default="combineResults", help="base output folder")
    parser.add_argument("-o", "--outfolder", type=str, default="", help="Main output folder, with the root file storing all histograms and datacards for single charge")
    parser.add_argument("-i", "--inputFile", type=str)
    parser.add_argument("--qcdScale", choices=["byHelicityPt", "byHelicityPtCharge", "byHelicityCharge", "byPtCharge", "byPt", "byCharge", "integrated",], default="byHelicityPtCharge", 
            help="Decorrelation for QCDscale")
    parser.add_argument("--rebinPtV", type=float, nargs='*', help="Rebin axis with gen boson pt by this value (default does nothing)")
    parser.add_argument("--scetlibUnc", action='store_true', help="Include SCETlib uncertainties")
    parser.add_argument("--pdf", type=str, default="nnpdf31", choices=theory_tools.pdfMapExtended.keys(), help="PDF to use")
    parser.add_argument("-b", "--fitObs", type=str, default="nominal", help="Observable to fit") # TODO: what does it do?
    parser.add_argument("--noQCDscaleFakes", dest="noQCDscaleFakes" , action="store_true",   help="Do not apply QCd scale uncertainties on fakes, mainly for debugging")
    parser.add_argument("--doStatOnly", action="store_true", default=False, help="Set up fit to get stat-only uncertainty (currently combinetf with -S 0 doesn't work)")
    parser.add_argument("--debug", action='store_true', help="Print debug output")
    return parser

def setup_base_logger(name, debug):
    logging.basicConfig()
    base_logger = logging.getLogger("wremnants")
    base_logger.setLevel(logging.DEBUG if debug else logging.INFO)
    return base_logger.getChild(name)
    
def child_logger(name):
    return logging.getLogger("wremnants").getChild(name.split(".")[-1])

'''
INPUT -------------------------------------------------------------------------
|* (str) string: the string to be converted to list
|
ROUTINE -----------------------------------------------------------------------
|* converts a string to a string element in a list
|  - if not comma-separated, then the whole string becomes one single element
OUTPUT ------------------------------------------------------------------------
|* (float) string: the list-lized string
+------------------------------------------------------------------------------
'''
def string_to_list(string):
	if type(string) == str:
		string = string.split(",") # items have to be comma-separated 
		return string
	elif type(string) == list:
		return string
	else:
		raise TypeError(
            "string_to_list(): cannot convert an input that is"
            "neither a single string nor a list of strings to a list"
        )

'''
INPUT -------------------------------------------------------------------------
|* list(str): a list of strings
|
ROUTINE -----------------------------------------------------------------------
|* convert the list of string to a single string by join()
|
OUTPUT ------------------------------------------------------------------------
|* (str): the resulted string
+------------------------------------------------------------------------------
'''
def list_to_string(list_str):
	if type(list_str) == str:
		return list_str
	elif type(list_str) == list:
		string = ""
		return string.join(list_str)
	else:
		raise TypeError(
            "list_to_string(): cannot convert an input that is"
            " neither a single string or a list of strings"
        )
