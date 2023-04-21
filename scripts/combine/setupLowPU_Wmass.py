#!/usr/bin/env python3
from wremnants import CardTool,theory_tools,syst_tools
from wremnants.datasets.datagroupsLowPU import make_datagroups_lowPU
from utilities import logging
from wremnants import histselections as sel
import argparse
import os
import pathlib
import hist
import numpy as np
import scripts.lowPU.config as lowPUcfg
import math

scriptdir = f"{pathlib.Path(__file__).parent}"

parser = argparse.ArgumentParser()
parser.add_argument("-o", "--outfolder", type=str, default="/scratch/jaeyserm/CombineStudies_Wmass")
parser.add_argument("-i", "--inputFile", type=str, default="")
parser.add_argument("--noScaleHelicitySplit", dest="qcdByHelicity", action='store_false', 
        help="Don't split QCD scale into helicity coefficients")
parser.add_argument("--qcdScale", choices=["byHelicityPt", "byPt", "integrated", "byHelicity"], default="integrated", 
        help="Decorrelation for QCDscale (additionally always by charge)")
parser.add_argument("--flavor", type=str, help="Flavor (e or mu)", default=None, required=True)
parser.add_argument("--fittype", choices=["differential", "wmass", "wlike", "inclusive"], default="differential", 
        help="Fit type, defines POI and fit observable (recoil or mT)")
parser.add_argument("--statOnly", action='store_true', help="Stat-only cards")
parser.add_argument("--lumiScale", help="Luminosity scale", type=float, default=1.0)
parser.add_argument("--met", type=str, help="MET (DeepMETReso or RawPFMET)", default="RawPFMET")
args = parser.parse_args()

logger = logging.setup_logger(__file__, args.verbose, args.noColorLogger)

if not os.path.isdir(args.outfolder):
    os.mkdir(args.outfolder)

met = args.met # mumu, ee
lumisuffix = str(args.lumiScale).replace(".", "p")
suffix = "_statOnly" if args.statOnly else ""

if args.inputFile == "": args.inputFile = "lowPU_%s_%s.pkl.lz4" % (args.flavor, met)

datagroups = make_datagroups_lowPU(args.inputFile, flavor=args.flavor)

unconstrainedProcs = ["WplusJetsToMuNu", "WminusJetsToMuNu"] # POIs
constrainedProcs = []   # constrained signal procs
bkgDataProcs = []       # all the rest + data
if args.flavor == "mu":
    dataName = "SingleMuon"
    bkgDataProcs = ["TTbar", "VV", "WJetsToTauNu", "DY", "Fake", "SingleMuon"]
if args.flavor == "ee":
    dataName = "SingleElectron"
    bkgDataProcs = ["TTbar", "EWK", "SingleElectron"]


histName = "mT_corr_rec"
QCDscale = "integral"


  
# hack: remove non-used procs/groups, as there can be more procs/groups defined than defined above
toDel = []
for group in datagroups.groups: 
    if not group in constrainedProcs+unconstrainedProcs+bkgDataProcs: toDel.append(group)
datagroups.deleteGroups(toDel)    

templateDir = f"{scriptdir}/Templates/LowPileupW"
cardTool = CardTool.CardTool(f"{args.outfolder}/LowPU_Wmass_{args.flavor}_{{chan}}_{met}_lumi{lumisuffix}{suffix}.txt")
cardTool.setNominalTemplate(f"{templateDir}/main.txt")
cardTool.setOutfile(os.path.abspath(f"{args.outfolder}/LowPU_Wmass_{args.flavor}_{met}_lumi{lumisuffix}{suffix}.root"))
cardTool.setDatagroups(datagroups)
cardTool.setHistName(histName) 
##cardTool.setChannels([f"{args.flavor}"])
cardTool.setDataName(dataName)
cardTool.setProcsNoStatUnc(procs=[])
cardTool.setSpacing(36)
##cardTool.setWriteByCharge(False)
cardTool.setLumiScale(args.lumiScale)
####cardTool.setUnconstrainedProcs(unconstrainedProcs)

logger.debug(f"Making datacards with these processes: {cardTool.getProcesses()}")

pdfName = theory_tools.pdfMapExtended["msht20"]["name"]
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

'''
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


recoil_vars = ["target_para", "target_perp", "source_para", "source_perp"]
#recoil_vars = ["target_para"]
for tag in recoil_vars:
    cardTool.addSystematic("recoilSyst_%s" % tag,
        processes=unconstrainedProcs,
        mirror = True,
        group = "CMS_recoil_stat",
        systAxes = ["axis_recoil_unc_%s" % tag],
        labelsByAxis = ["recoilSyst_%s_{i}" % tag],
        #action=scale_recoil_hist_to_variations,
        #scale = 1./math.sqrt(args.lumiScale),
    )


# dummy muon momentum scale
cardTool.addSystematic("muonScaleSyst", 
    processes=unconstrainedProcs,
    group="muonScale",
    baseName="CMS_scale_m_",
    systAxes=["downUpVar", "scaleEtaSlice"],
    labelsByAxis=["downUpVar", "ieta"],
    scale = 0.5, # half of the uncertainty propagates to mT
)

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
    

cardTool.addLnNSystematic("CMS_Top", processes=["TTbar"], size=1.06, group="CMS_Top")
cardTool.addLnNSystematic("CMS_VV", processes=["VV"], size=1.16, group="CMS_VV")
cardTool.addLnNSystematic("CMS_DY", processes=["DY"], size=1.03, group="CMS_DY")
#cardTool.addLnNSystematic("CMS_WJetsToTauNu", processes=["WplusJetsToTauNu", "WminusJetsToTauNu"], size=1.03, group="CMS_WJetsToTauNu")
#cardTool.addLnNSystematic("CMS_WjetsMuNu", processes=unconstrainedProcs, size=1.06, group="CMS_WjetsMuNu")
cardTool.addLnNSystematic("CMS_lumi_lowPU", processes=cardTool.allMCProcesses(), size=1.02, group="CMS_lumi_lowPU")


cardTool.writeOutput(args=args)
