#!/usr/bin/env python3
from wremnants import CardTool,theory_tools,syst_tools

from wremnants import histselections as sel
import argparse
import os
import pathlib
import hist
import numpy as np
import scripts.lowPU.config as lowPUcfg
import math
from utilities import common, logging

scriptdir = f"{pathlib.Path(__file__).parent}"

parser = argparse.ArgumentParser()
parser.add_argument("-o", "--outfolder", type=str, default="/scratch/jaeyserm/CombineStudies_Wmass_mT")
parser.add_argument("--noScaleHelicitySplit", dest="qcdByHelicity", action='store_false', 
        help="Don't split QCD scale into helicity coefficients")
#parser.add_argument("--qcdScale", choices=["byHelicityPt", "byPt", "integrated", "byHelicity"], default="integrated", 
#        help="Decorrelation for QCDscale (additionally always by charge)")
parser.add_argument("--flavor", type=str, help="Flavor (e or mu)", default="mu")
parser.add_argument("--fittype", choices=["differential", "wmass", "wlike", "inclusive"], default="differential", 
        help="Fit type, defines POI and fit observable (recoil or mT)")
parser.add_argument("--statOnly", action='store_true', help="Stat-only cards")
parser.add_argument("--lumiScale", help="Luminosity scale", type=float, default=1.0)
parser.add_argument("--met", type=str, help="MET (DeepMETReso or RawPFMET)", default="RawPFMET")
parser.add_argument("--PUMode", type=str, help="PU mode (lowPU or highPU)", default="lowPU")

parser.add_argument("--qcdScale", choices=["byHelicityPtAndByPt", "byHelicityPt", "byHelicityCharge", "byPt", "byCharge", "integrated",], default="integrated", 
        help="Decorrelation for QCDscale (additionally always by charge). With 'byHelicityPtAndByPt' two independent histograms are stored, split and not split by helicities (for tests)")
parser.add_argument("--noQCDscaleFakes", action="store_true",   help="Do not apply QCd scale uncertainties on fakes, mainly for debugging")
parser.add_argument("--skipSignalSystOnFakes", action="store_true", help="Do not propagate signal uncertainties on fakes, mainly for checks.")


args = parser.parse_args()
logger = logging.setup_logger(__file__, args.verbose, args.noColorLogger)

if not os.path.isdir(args.outfolder):
    os.mkdir(args.outfolder)

met = args.met # mumu, ee
lumisuffix = str(args.lumiScale).replace(".", "p")
suffix = "_statOnly" if args.statOnly else ""

if args.PUMode == "lowPU":

    from wremnants.datasets.datagroupsLowPU import make_datagroups_lowPU
    inputFile = "lowPU_%s_%s.pkl.lz4" % (args.flavor, met)
    datagroups = make_datagroups_lowPU(inputFile, flavor=args.flavor)
    
    #unconstrainedProcs = ["WplusJetsToMuNu", "WminusJetsToMuNu"] # POIs
    unconstrainedProcs = ["WJetsToMuNu"] # POIs
    
    constrainedProcs = []   # constrained signal procs
    bkgDataProcs = []       # all the rest + data
    dataName = "SingleMuon"
    bkgDataProcs = ["TTbar", "VV", "WJetsToTauNu", "DY", "Fake", "SingleMuon"]

elif args.PUMode == "highPU":

    from wremnants.datasets.datagroups import datagroups2016_mT
    inputFile = "mw_with_%s_eta_pt_%s.pkl.lz4" % (args.flavor, met)
    datagroups = datagroups2016_mT(inputFile)
    
    unconstrainedProcs = ["Wmunu"] # POIs
    constrainedProcs = []   # constrained signal procs
    bkgDataProcs = []       # all the rest + data
    dataName = "Data"
    bkgDataProcs = ["Top", "EWK",  "Data"]
    
    ## python scripts/combine/setupCombineWMass_mT.py --PUMode=highPU --met=RawPFMET --statOnly
    ## python scripts/combine/setupCombineWMass_mT.py --PUMode=highPU --met=DeepMETReso --statOnly

else: 
    print("PU mode should be either lowPU or highPU")
    quit()




histName = "mT_corr_rec"
QCDscale = "integral"


  
# hack: remove non-used procs/groups, as there can be more procs/groups defined than defined above
toDel = []
for group in datagroups.groups: 
    if not group in constrainedProcs+unconstrainedProcs+bkgDataProcs: toDel.append(group)
datagroups.deleteGroups(toDel)    

templateDir = f"{scriptdir}/Templates/LowPileupW"
cardTool = CardTool.CardTool(f"{args.outfolder}/{args.PUMode}_{args.flavor}_{{chan}}_{met}_lumi{lumisuffix}{suffix}.txt")
cardTool.setNominalTemplate(f"{templateDir}/main.txt")
cardTool.setOutfile(os.path.abspath(f"{args.outfolder}/{args.PUMode}_{args.flavor}_{met}_lumi{lumisuffix}{suffix}.root"))
cardTool.setDatagroups(datagroups)
cardTool.setHistName(histName) 
##cardTool.setChannels([f"{args.flavor}"])
cardTool.setDataName(dataName)
cardTool.setProcsNoStatUnc(procs=[])
cardTool.setSpacing(36)
##cardTool.setWriteByCharge(False)
cardTool.setLumiScale(args.lumiScale)
####cardTool.setUnconstrainedProcs(unconstrainedProcs)


if args.statOnly:
    # print a card with only mass weights and a dummy syst
    cardTool.addLnNSystematic("dummy", processes=cardTool.allMCProcesses(), size=1.00001, group="dummy")
    
    cardTool.addSystematic("massWeight", 
        processes=unconstrainedProcs,
        outNames=theory_tools.massWeightNames(["massShift100MeV"], wlike=False),
        group="massShift",
        groupFilter=lambda x: x == "massShift100MeV",
        mirror=False,
        #TODO: Name this
        noConstraint=True,
        systAxes=["tensor_axis_0"],
    )
    
    
    cardTool.writeOutput(args=args)
    print("Using option --doStatOnly: the card was created with only mass weights and a dummy LnN syst on all processes")
    quit()

'''
pdfName = theory_tools.pdfMap["nnpdf31"]["name"]
cardTool.addSystematic(pdfName, 
    processes=unconstrainedProcs,
    mirror=True,
    group=pdfName,
    systAxes=["tensor_axis_0"],
    labelsByAxis=[pdfName.replace("pdf", "pdf{i}")],
    # Needs to be a tuple, since for multiple axis it would be (ax1, ax2, ax3)...
    # -1 means all possible values of the mirror axis
    skipEntries=[(0, -1)],
)
cardTool.addSystematic(f"alphaS002{pdfName}", 
    processes=unconstrainedProcs,
    mirror=False,
    group=pdfName,
    systAxes=["tensor_axis_0"],
    outNames=[pdfName+"AlphaSUp", pdfName+"AlphaSDown"],
    scale=0.75,
)


### QCD SCALES


scaleSystAxes = ["muRfact", "muFfact"] 
scaleLabelsByAxis = ["muR", "muF"]
scaleGroupName = "QCDscale"
# Exclude all combinations where muR = muF = 1 (nominal) or where
# they are extreme values (ratio = 4 or 1/4)
scaleSkipEntries = [(1, 1), (0, 2), (2, 0)]
# This is hacky but it's the best idea I have for now...
systNameReplaceVec = [("muR2muF2", "muRmuFUp"), ("muR0muF0", "muRmuFDown"), ("muR2muF1", "muRUp"), 
                      ("muR0muF1", "muRDown"), ("muR1muF0", "muFDown"), ("muR1muF2", "muFUp"),
                      ("genQ0", "genVminus"), ("genQ1", "genVplus")]


inclusiveScale = args.qcdScale == "integrated"
helicity = "Helicity" in args.qcdScale
passSystToFakes = False # not args.wlike and not args.skipSignalSystOnFakes

if inclusiveScale:
    scale_hist = "inclusive_qcdScale"
    scale_action = None
    scaleActionArgs = {}
    scale_action_map = {}
else:
    scale_action_map = {proc : syst_tools.scale_helicity_hist_to_variations for proc in (common.zprocs if args.wlike else common.wprocs)}
    scale_hist = "qcdScale"
    
if args.qcdScale == "byCharge":
    scaleActionArgs = {"sum_axis" : ["ptVgen"]} 
    scaleGroupName += "ByChargeV"
    scaleSystAxes.insert(0, "chargeVgen")
    scaleLabelsByAxis.insert(0, "genQ")
    scaleSkipEntries = [(-1, *x) for x in scaleSkipEntries] # need to add a -1 for each axis element added before
    
if "Pt" in args.qcdScale:
    scaleActionArgs = {"rebinPtV" : args.rebinPtV}
    scaleGroupName += "ByPtV"
    scaleSystAxes.insert(0, "ptVgen")
    scaleLabelsByAxis.insert(0, "PtBin")
    scaleSystAxes.insert(0, "chargeVgen")
    scaleLabelsByAxis.insert(0, "genQ")
    scaleSkipEntries = [(-1, -1, *x) for x in scaleSkipEntries] # need to add a -1 for each axis element added before

    ## keep commented for now, it is only for tests but might become useful to deal with constraints
    ##
    # if args.qcdScale == "byHelicityPtAndByPt":
    #     # note: the following uses a different histogram compared to the one for helicity splitting
    #     # when integrating on the coefficients for the helicity split version they should give the same results modulo statistical fluctuations
    #     print("Option --qcdScale byHelicityPtAndByPt was chosen: doing additional qcdScale histogram for tests")
    #     print(scaleActionArgs)
    #     print(scaleLabelsByAxis)
    #     cardTool.addSystematic("qcdScale",
    #                            actionMap=scale_action_map,
    #                            actionArgs=copy.deepcopy(scaleActionArgs), # to avoid possible undesired updates below
    #                            processes=signal_samples_inctau,
    #                            group=scaleGroupName,
    #                            systAxes=scaleSystAxes[:],
    #                            labelsByAxis=scaleLabelsByAxis[:],
    #                            skipEntries=scaleSkipEntries[:],
    #                            systNameReplace=systNameReplaceVec,
    #                            baseName="QCDscaleByPt_",
    #                            passToFakes=False if args.noQCDscaleFakes else passSystToFakes,
    #     )

if helicity:
    scale_hist = "qcdScaleByHelicity"
    scaleActionArgs = {"rebinPtV" : args.rebinPtV}
    scaleSystAxes.insert(0, "helicity")
    scaleLabelsByAxis.insert(0, "Coeff")
    if args.qcdScale == "byHelicityCharge":
        # mainly for tests
        scaleActionArgs.update({"sum_axis" : ["ptVgen"]})
        scaleGroupName += "ByHelicityCharge"
        scaleSystAxes.insert(0, "chargeVgen")
        scaleLabelsByAxis.insert(0, "genQ")
        scaleSkipEntries = [(-1, -1, *x) for x in scaleSkipEntries] # need to add a -1 for each axis element added before
    else:
        scaleGroupName += "ByHelicity"
        scaleSkipEntries = [(-1, *x) for x in scaleSkipEntries] # need to add a -1 for each axis element added before

print("Inclusive scale", inclusiveScale)
print(scaleActionArgs if not inclusiveScale else None)
print(scaleLabelsByAxis)

scaleSkipEntries = [(-1, 1, 1), (-1, 0, 2), (-1, 2, 0)] 
scaleSystAxes = ["chargeVgen", "muRfact", "muFfact"]
scaleLabelsByAxis = ["q", "muR", "muF"]
scale_action_map = {proc : syst_tools.scale_helicity_hist_to_variations for proc in ["WplusJetsToMuNu", "WminusJetsToMuNu", "WplusJetsToTauNu", "WminusJetsToTauNu"]} # must be the proc name, not group name

cardTool.addSystematic("qcdScaleByHelicity",
    actionMap=scale_action_map,
    actionArgs={"sum_axis" : ["ptVgen", "helicity"]},
    ###processes=signal_samples_inctau,
    processes=unconstrainedProcs + ["WJetsToTauNu"], #+ ["WJetsToTauNu"],
    group=scaleGroupName,
    # splitGroup={f"{scaleGroupName}_coeff{i}" : f".*Coeff{i}" for i in range(9)}, # key is the new group name to make it unique, value is the pattern to filter nuisances
    systAxes=scaleSystAxes,
    labelsByAxis=scaleLabelsByAxis,
    # Exclude all combinations where muR = muF = 1 (nominal) or where
    # they are extreme values (ratio = 4 or 1/4)
    skipEntries=scaleSkipEntries,
    systNameReplace=systNameReplaceVec,
    baseName="QCDscale_",
    passToFakes=False,
)

'''

'''
scaleSystAxes = ["chargeVgen", "muRfact", "muFfact"]  ##  "ptVgen", "chargeVgen", "helicityWeight_tensor" ## charge = 0?
scaleLabelsByAxis = ["q", "muR", "muF"]
scaleGroupName = "QCDscale"
scaleActionArgs = {"sum_ptV" : True, "sum_helicity" : True}
scaleSkipEntries = [(-1, 1, 1), (-1, 0, 2), (-1, 2, 0)] 
if "Helicity" in QCDscale:
    scaleActionArgs.pop("sum_helicity")
    scaleGroupName += "ByHelicity"
    scaleSystAxes.insert(0, "helicity")
    scaleLabelsByAxis.insert(0, "Coeff")
    scaleSkipEntries = [(-1, *x) for x in scaleSkipEntries]
if "Pt" in QCDscale:
    scaleActionArgs.pop("sum_ptV")
    scaleGroupName += "ByPtV"
    scaleSystAxes.insert(0, "ptVgen")
    scaleLabelsByAxis.insert(0, "genPtV")
    scaleSkipEntries = [(-1, *x) for x in scaleSkipEntries]

cardTool.addSystematic("qcdScaleByHelicity", 
    action=syst_tools.scale_helicity_hist_to_variations,
    actionArgs=scaleActionArgs,
    processes=unconstrainedProcs + ["WJetsToTauNu"],
    group=scaleGroupName,
    systAxes=scaleSystAxes,
    labelsByAxis=scaleLabelsByAxis,
    # Exclude all combinations where muR = muF = 1 (nominal) or where
    # they are extreme values (ratio = 4 or 1/4)
    skipEntries=scaleSkipEntries,
    # This is hacky but it's the best idea I have for now...
    systNameReplace=[("muR2muF2", "muRmuFUp"), ("muR0muF0", "muRmuFDown"), ("muR2muF1", "muRUp"), 
        ("muR0muF1", "muRDown"), ("muR1muF0", "muFDown"), ("muR1muF2", "muFUp")],
    baseName="QCDscale_",
    )
'''

'''
def scale_recoil_hist_to_variations(scale_hist):
    ## "mt_recoilStatUnc", "mt", "mt_pert"
    # scale_hist = recoil_gen, recoil_reco, recoil_reco_pert, mll, systIdx (=qT bin variations)
    # recoil_gen already removed
    
    #if hName == "recoil_corr_rec_para_qT": basname_syst, axis, axis_pert = "recoil_para_qT_recoilStatUnc", "recoil_para_qT", "recoil_para_qT_pert"
    #elif hName == "recoil_corr_rec_perp": basname_syst, axis, axis_pert = "recoil_perp_recoilStatUnc", "recoil_perp", "recoil_perp_pert"
    #elif hName == "MET_corr_rec_pt": basname_syst, axis, axis_pert = "MET_recoilStatUnc", "recoil_MET_pt", "recoil_MET_pt_pert"
    #elif hName == "mT_corr_rec": basname_syst, axis, axis_pert = "mt_recoilStatUnc", "mt", "mt_pert"
    #elif hName == "recoil_corr_rec_magn": basname_syst, axis, axis_pert = "recoil_magn_recoilStatUnc", "recoil_magn", "recoil_magn_pert"
    basname_syst, axis, axis_pert = "mt_recoilStatUnc", "mt", "mt_pert"

    # BASE: sum over qT bins, RECO only --> (RECO, M)
    # NOM: RECO only, in bins of qT --> (RECO, M, QT)
    # PERT: RECOpert only, in bins of qT --> (RECOpert, M, QT)
    # --> VAR = NOM(extended) - NOM + PERT
    
    s = hist.tag.Slicer()
    #base = scale_hist[{"recoil_stat_unc_var" : s[::hist.sum], "recoil_reco_pert" : s[::hist.sum]}] # sum over qT, remove the qT axis == (RECO, M)
    #nom = scale_hist[{"recoil_reco_pert" : s[::hist.sum]}] # sum over qT (perturbed)  == (RECO, M, IDX)
    #pert = scale_hist[{"recoil_reco" : s[::hist.sum]}] # sum over qT == (RECO_PERT, M, IDX)
    base = scale_hist[{"recoil_stat_unc_var" : s[::hist.sum], axis_pert : s[::hist.sum]}] # sum over qT, remove the qT axis == (RECO, M)
    nom = scale_hist[{axis_pert : s[::hist.sum]}] # sum over qT (perturbed)  == (RECO, M, IDX)
    pert = scale_hist[{axis : s[::hist.sum]}] # sum over qT == (RECO_PERT, M, IDX)
    
    #print(base)
    #print(nom)
    #print(pert)
    #print(scale_hist)
    
    out_name = "recoilStat_variations"
    
    # up
    # 4 D hist
    #  base.view(flow=True)[..., np.newaxis] --> indentical copies of the nominal
    scale_variation_hist_up = hist.Hist(*pert.axes, storage = scale_hist._storage_type(), name = out_name+"Up",
                data = pert.view(flow=True) - nom.view(flow=True) + base.view(flow=True)[..., np.newaxis, :, :])
       
    scale_variation_hist_dw = hist.Hist(*pert.axes, storage = scale_hist._storage_type(), name = out_name+"Down",
                data = - scale_variation_hist_up.view(flow=True) +  base.view(flow=True)[..., np.newaxis, :, :] +  base.view(flow=True)[..., np.newaxis, :, :] )
       
    

                 
    #hmirror = hh.mirrorHist(scale_variation_hist_up, nom)
    mirrorAx = hist.axis.Integer(0,2, name="mirror", overflow=False, underflow=False)
    hnew = hist.Hist(*scale_variation_hist_up.axes, mirrorAx, storage=scale_variation_hist_up._storage_type(), name = out_name)
    hnew.view(flow=True)[...] = np.stack((scale_variation_hist_up.view(flow=True), scale_variation_hist_dw.view(flow=True)), axis=-1)
    
    return hnew      


recoil_vars = [(1,2), (1,3), (1,4),   (2,2), (2,3),   (3,2), (3,3), (3,4),    (4,2), (4,3)]
recoil_vars = [(1,2),(1,3),(1,4),(1,5),  (2,2),(2,3),(2,4),  (3,2),(3,3),(3,4),(3,5),  (4,2),(4,3),(4,4)]
#recoil_vars = [(1,2)]
#recoil_vars = []
for k in recoil_vars:
    
    cardTool.addSystematic("recoilStatUnc_%d_%d" % (k[0], k[1]),
        processes=unconstrainedProcs,
        mirror = False, # mirroring done in scale_recoil_hist_to_variations()
        group = "CMS_recoil_stat",
        systAxes = ["recoil_stat_unc_var"],
        labelsByAxis = ["recoilStatUnc_%d_%d_{i}" % (k[0], k[1])],
        action=scale_recoil_hist_to_variations,
        scale = 1./math.sqrt(args.lumiScale),
    )

'''

'''
recoil_vars = ["target_para", "target_perp", "source_para", "source_perp", "target_para_bkg", "target_perp_bkg"]
recoil_grps = ["recoil_stat", "recoil_stat", "recoil_stat", "recoil_stat", "recoil_syst", "recoil_syst"]
for i, tag in enumerate(recoil_vars):
    cardTool.addSystematic("recoilSyst_%s" % tag,
        processes=unconstrainedProcs,
        mirror = True,
        group = recoil_grps[i],
        systAxes = ["axis_recoil_unc_%s" % tag],
        labelsByAxis = ["recoil_%s_{i}" % tag],
        #action=scale_recoil_hist_to_variations,
        #scale = 1./math.sqrt(args.lumiScale),
    )


'''
recoil_vars = ["target_para", "target_perp", "source_para", "source_perp", "target_para_bkg", "target_perp_bkg"]
recoil_grps = ["recoil_stat", "recoil_stat", "recoil_stat", "recoil_stat", "recoil_syst", "recoil_syst"]
for i, tag in enumerate(recoil_vars):
    cardTool.addSystematic("recoilSystWeight_%s" % tag,
        processes=unconstrainedProcs,
        mirror = True,
        group = recoil_grps[i],
        systAxes = ["tensor_axis_0"],
        labelsByAxis = ["recoil_%s_{i}" % tag],
        #action=scale_recoil_hist_to_variations,
        #scale = 1./math.sqrt(args.lumiScale) if not "bkg" in tag else 1.0,
    )

'''
# dummy muon momentum scale
cardTool.addSystematic("muonScaleSyst", 
    processes=unconstrainedProcs,
    group="muonScale",
    baseName="CMS_scale_m_",
    systAxes=["downUpVar", "scaleEtaSlice"],
    labelsByAxis=["downUpVar", "ieta"],
    scale = 0.5, # half of the uncertainty propagates to mT
)
'''



cardTool.addSystematic("massWeight", 
    processes=unconstrainedProcs,
    outNames=theory_tools.massWeightNames(["massShift100MeV"], wlike=False),
    group="massShift",
    groupFilter=lambda x: x == "massShift100MeV",
    mirror=False,
    #TODO: Name this
    noConstraint=True,
    systAxes=["tensor_axis_0"],
)
    

#cardTool.addLnNSystematic("dummy", processes=cardTool.allMCProcesses(), size=1.001, group="dummy")
###cardTool.addLnNSystematic("CMS_Top", processes=["TTbar"], size=1.06, group="CMS_Top")
###cardTool.addLnNSystematic("CMS_VV", processes=["VV"], size=1.16, group="CMS_VV")
###cardTool.addLnNSystematic("CMS_DY", processes=["DY"], size=1.03, group="CMS_DY")
###cardTool.addLnNSystematic("CMS_Fakes", processes=["Fake"], size=1.05, group="CMS_Fakes_norm")
#cardTool.addLnNSystematic("CMS_WJetsToTauNu", processes=["WplusJetsToTauNu", "WminusJetsToTauNu"], size=1.03, group="CMS_WJetsToTauNu")
#cardTool.addLnNSystematic("CMS_WjetsMuNu", processes=unconstrainedProcs, size=1.06, group="CMS_WjetsMuNu")
###cardTool.addLnNSystematic("CMS_lumi_lowPU", processes=cardTool.allMCProcesses(), size=1.02, group="CMS_lumi_lowPU")


cardTool.writeOutput(args=args)
