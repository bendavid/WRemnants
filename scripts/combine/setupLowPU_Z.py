#!/usr/bin/env python3
from wremnants import CardTool,theory_tools,syst_tools
from wremnants.datasets.datagroupsLowPU import datagroupsLowPU_Z
from wremnants import histselections as sel
import argparse
import os
import pathlib
import hist
import numpy as np
import scripts.lowPU.config as lowPUcfg

scriptdir = f"{pathlib.Path(__file__).parent}"

parser = argparse.ArgumentParser()
parser.add_argument("-o", "--outfolder", type=str, default="/scratch/jaeyserm/CombineStudies")
parser.add_argument("-i", "--inputFile", type=str, default="")
parser.add_argument("--noScaleHelicitySplit", dest="qcdByHelicity", action='store_false', 
        help="Don't split QCD scale into helicity coefficients")
parser.add_argument("--qcdScale", choices=["byHelicityPt", "byPt", "integrated", "byHelicity"], default="integrated", 
        help="Decorrelation for QCDscale (additionally always by charge)")
parser.add_argument("--flavor", type=str, help="Flavor (ee or mumu)", default=None, required=True)
parser.add_argument("--fittype", choices=["differential", "wmass", "wlike", "inclusive"], default="differential", 
        help="Fit type, defines POI and fit observable (recoil or mT)")
parser.add_argument("--xsec", dest="xsec", action='store_true', 
        help="Write card for masked xsec normalization")
args = parser.parse_args()

if not os.path.isdir(args.outfolder):
    os.mkdir(args.outfolder)

if args.inputFile == "": args.inputFile = "mz_lowPU_%s.pkl.lz4" % args.flavor

datagroups = datagroupsLowPU_Z(args.inputFile, flavor=args.flavor)

unconstrainedProcs = [] # POIs
constrainedProcs = []   # constrained signal procs
bkgDataProcs = []       # all the rest + data
suffix = ""
if args.flavor == "mumu":
    dataName = "SingleMuon"
    bkgDataProcs = ["TTbar", "EWK", "SingleMuon"] # , "DYmumu"
if args.flavor == "ee":
    dataName = "SingleElectron"
    bkgDataProcs = ["TTbar", "EWK", "SingleElectron"] # , "DYee"
if args.xsec:
    suffix = "_xsec"
    bkgDataProcs = [bkgDataProcs[-1]] # keep data_obs

recoBins = lowPUcfg.bins_recoil_reco
genBins = lowPUcfg.bins_recoil_gen


histName = ""
QCDscale = ""
if args.fittype == "differential" or args.fittype == "wmass":
    histName = "reco_mll" # for signal gen_reco_mll

    proc_base = dict(datagroups.groups["DYmumu" if args.flavor == "mumu" else "DYee"]) # load the base process (DYmumu or DYee)
    for i in range(len(genBins)-1): # add gen bin processes
    
        proc_name = "DY_genBin%d" % (i+1)
        proc_genbin = dict(proc_base)
        proc_genbin['signalOp'] = lambda x, i=i: x[{"recoil_gen" : i}] # index correct? Start from 1?
        datagroups.groups[proc_name] = proc_genbin
        if args.fittype == "differential": unconstrainedProcs.append(proc_name)
        elif args.fittype == "wmass": constrainedProcs.append(proc_name)
    if args.fittype == "differential": QCDscale = "integral"
    elif args.fittype == "wmass": QCDscale = "byHelicityPt"

elif args.fittype == "wlike":
    histName = "mt"
    constrainedProcs.append("DYmumu" if args.flavor == "mumu" else "DYee")
    QCDscale = "integral"
elif args.fittype == "inclusive":
    histName = "reco_mll"
    unconstrainedProcs.append("DYmumu" if args.flavor == "mumu" else "DYee")
    QCDscale = "integral"
    

  
# hack: remove non-used procs/groups, as there can be more procs/groups defined than defined above
toDel = []
for group in datagroups.groups: 
    if not group in constrainedProcs+unconstrainedProcs+bkgDataProcs: toDel.append(group)
for group in toDel: del datagroups.groups[group]    

templateDir = f"{scriptdir}/Templates/LowPileupW"
cardTool = CardTool.CardTool(f"{args.outfolder}/LowPU_Z{args.flavor}_{args.fittype}{suffix}.txt")
cardTool.setNominalTemplate(f"{templateDir}/main.txt")
cardTool.setOutfile(os.path.abspath(f"{args.outfolder}/LowPU_Z{args.flavor}_{args.fittype}{suffix}.root"))
cardTool.setDatagroups(datagroups)
cardTool.setHistName(histName) 
cardTool.setChannels([f"{args.flavor}{suffix}"])
cardTool.setDataName(dataName)
cardTool.setUnconstrainedProcs(unconstrainedProcs)

DY_procs = cardTool.filteredProcesses(lambda x: "DY" in x)


pdfName = theory_tools.pdfMap["nnpdf31"]["name"]
cardTool.addSystematic(pdfName, 
    processes=DY_procs,
    mirror=True,
    group=pdfName,
    systAxes=["tensor_axis_0"],
    labelsByAxis=[pdfName.replace("pdf", "pdf{i}")],
    # Needs to be a tuple, since for multiple axis it would be (ax1, ax2, ax3)...
    # -1 means all possible values of the mirror axis
    skipEntries=[(0, -1)],
)
cardTool.addSystematic(f"alphaS002{pdfName}", 
    processes=DY_procs,
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
    processes=DY_procs,
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

    # scale_hist = recoil_gen, recoil_reco, recoil_reco_pert, mll, systIdx (=qT bin variations)
    # recoil_gen already removed

    # BASE: sum over qT bins, RECO only --> (RECO, M)
    # NOM: RECO only, in bins of qT --> (RECO, M, QT)
    # PERT: RECOpert only, in bins of qT --> (RECOpert, M, QT)
    # --> VAR = NOM(extended) - NOM + PERT
    
    s = hist.tag.Slicer()
    base = scale_hist[{"recoil_stat_unc_var" : s[::hist.sum], "recoil_reco_pert" : s[::hist.sum]}] # sum over qT, remove the qT axis == (RECO, M)
    nom = scale_hist[{"recoil_reco_pert" : s[::hist.sum]}] # sum over qT (perturbed)  == (RECO, M, IDX)
    pert = scale_hist[{"recoil_reco" : s[::hist.sum]}] # sum over qT == (RECO_PERT, M, IDX)
    
    #print(scale_hist)
    
    out_name = "recoilStat_variations"
    
    # up
    # 4 D hist
    #  base.view(flow=True)[..., np.newaxis] --> indentical copies of the nominal
    scale_variation_hist_up = hist.Hist(*pert.axes, storage = scale_hist._storage_type(), name = out_name+"Up",
                data = pert.view(flow=True) - nom.view(flow=True) + base.view(flow=True)[..., np.newaxis])
       
    scale_variation_hist_dw = hist.Hist(*pert.axes, storage = scale_hist._storage_type(), name = out_name+"Down",
                data = - scale_variation_hist_up.view(flow=True) +  base.view(flow=True)[..., np.newaxis] +  base.view(flow=True)[..., np.newaxis] )
       
    

                 
    #hmirror = hh.mirrorHist(scale_variation_hist_up, nom)
    mirrorAx = hist.axis.Integer(0,2, name="mirror", overflow=False, underflow=False)
    hnew = hist.Hist(*scale_variation_hist_up.axes, mirrorAx, storage=scale_variation_hist_up._storage_type(), name = out_name)
    hnew.view(flow=True)[...] = np.stack((scale_variation_hist_up.view(flow=True), scale_variation_hist_dw.view(flow=True)), axis=-1)
    
    return hnew      


if not args.xsec:
   
    recoil_vars = [(1,2), (1,3), (1,4),   (2,2), (2,3),   (3,2), (3,3), (3,4),    (4,2), (4,3)]
    for k in recoil_vars:
        
        cardTool.addSystematic("recoilStatUnc_%d_%d" % (k[0], k[1]),
            processes=DY_procs,
            mirror = False,
            group = "CMS_recoil_stat",
            systAxes = ["recoil_stat_unc_var"],
            labelsByAxis = ["recoilStatUnc_%d_%d_{i}" % (k[0], k[1])],
            action=scale_recoil_hist_to_variations,
        )


    cardTool.addSystematic("prefireCorr",
        processes=DY_procs,
        mirror = False,
        group="CMS_prefire17",
        baseName="CMS_prefire17",
        systAxes = ["downUpVar"],
        labelsByAxis = ["downUpVar"],
    )
    
    for lepEff in ["lepSF_HLT_DATA_stat", "lepSF_HLT_DATA_syst", "lepSF_HLT_MC_stat", "lepSF_HLT_MC_syst", "lepSF_ISO_stat", "lepSF_ISO_DATA_syst", "lepSF_ISO_MC_syst", "lepSF_IDIP_stat", "lepSF_IDIP_DATA_syst", "lepSF_IDIP_MC_syst"]:
        
        
        cardTool.addSystematic(lepEff,
            processes=DY_procs,
            mirror = True,
            group="CMS_lepton_eff",
            baseName=lepEff,
            systAxes = ["tensor_axis_0"],
            labelsByAxis = [""],
        )
        
    cardTool.addLnNSystematic("CMS_Top", processes=["TTbar"], size=1.06, group="CMS_bkg_norm")
    cardTool.addLnNSystematic("CMS_VV", processes=["EWK"], size=1.16, group="CMS_bkg_norm")
    cardTool.addLnNSystematic("CMS_lumi_lowPU", processes=cardTool.allMCProcesses(), size=1.02, group="CMS_lumi_lowPU")


if args.fittype == "wlike" or args.fittype == "wmass":

    cardTool.addSystematic("massWeight", 
        processes=constrainedProcs,
        outNames=theory_tools.massWeightNames(["massShift100MeV"], wlike=False),
        group="massShift",
        groupFilter=lambda x: x == "massShift100MeV",
        mirror=False,
        #TODO: Name this
        noConstraint=True,
        systAxes=["tensor_axis_0"],
    )
    


cardTool.writeOutput()
