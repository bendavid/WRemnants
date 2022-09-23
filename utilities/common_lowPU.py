import hist
import pathlib
import argparse
import logging
import math

wremnants_dir = f"{pathlib.Path(__file__).parent}/../wremnants"
data_dir = f"{wremnants_dir}/data/"

wprocs = ["WminusJetsToMuNu", "WminusJetsToENu", "WminusJetsToTauNu", "WplusJetsToMuNu", "WplusJetsToENu", "WplusJetsToTauNu"]
zprocs = ["Zmumu", "Zee", "Ztautau"]
zprocs_recoil = ["Zmumu", "Zee"]
vprocs = wprocs+zprocs


# unfolding axes
axis_recoil_reco_ptZ = hist.axis.Variable([0, 5, 10, 15, 20, 30, 40, 50, 60, 75, 90, 150], name = "recoil_reco", underflow=False, overflow=True)
axis_recoil_gen_ptZ = hist.axis.Variable([0.0, 10.0, 20.0, 40.0, 60.0, 90.0, 150], name = "recoil_gen", underflow=False, overflow=True)


# categorical axes in python bindings always have an overflow bin, so use a regular
# axis for the charge
axis_charge = hist.axis.Regular(2, -2., 2., underflow=False, overflow=False, name = "charge")
axis_passIso = hist.axis.Boolean(name = "passIso")
axis_passMT = hist.axis.Boolean(name = "passMT")
down_up_axis = hist.axis.Regular(2, -2., 2., underflow=False, overflow=False, name = "downUpVar")






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
    
    parser.add_argument("--flavor", type=str, help="Flavor (ee or mumu)", default=None)
    parser.add_argument("--met", type=str, help="MET (DeepMETReso or RawPFMET)", default="RawPFMET")

    parser.add_argument("--pdfs", type=str, nargs="*", default=["nnpdf31"], choices=theory_tools.pdfMapExtended.keys(), help="PDF sets to produce error hists for")
    parser.add_argument("--altPdfOnlyCentral", action='store_true', help="Only store central value for alternate PDF sets")
    parser.add_argument("--maxFiles", type=int, help="Max number of files (per dataset)", default=-1)
    parser.add_argument("--filterProcs", type=str, nargs="*", help="Only run over processes matched by (subset) of name", default=[])
    parser.add_argument("--v8", action='store_true', help="Use NanoAODv8. Default is v9")
    parser.add_argument("-p", "--postfix", type=str, help="Postfix for output file name", default=None)
    parser.add_argument("--theory_corr", nargs="*", choices=["scetlib", "scetlibMSHT20", "scetlibHelicity", "dyturbo", "matrix_radish"], 
        help="Apply corrections from indicated generator. First will be nominal correction.", default=[])
    parser.add_argument("--theory_corr_alt_only", action='store_true', help="Save hist for correction hists but don't modify central weight")
    parser.add_argument("--skipHelicity", action='store_true', help="Skip the qcdScaleByHelicity histogram (it can be huge)")
    parser.add_argument("--eta", nargs=3, type=float, help="Eta binning as 'nbins min max' (only uniform for now)", default=[48,-2.4,2.4])
    parser.add_argument("--pt", nargs=3, type=float, help="Pt binning as 'nbins,min,max' (only uniform for now)", default=[29,26.,55.])
    parser.add_argument("--no_recoil", action='store_true', help="Don't apply recoild correction")
    return parser,initargs
