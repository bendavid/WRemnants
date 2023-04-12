#!/usr/bin/env python3
from wremnants import CardTool,theory_tools,syst_tools,combine_helpers
from wremnants.datasets.datagroupsLowPU import datagroupsLowPU
from utilities import common, logging
import argparse
import os
import pathlib
import hist
import numpy as np

scriptdir = f"{pathlib.Path(__file__).parent}"

def make_parser(parser=None):
    if not parser:
        parser = common.common_parser_combine()
    parser.add_argument("--flavor", choices=["ee", "mumu"], help="Flavor (ee or mumu)", default="mumu")
    parser.add_argument("--fitType", choices=["differential", "wmass", "wlike", "inclusive"], default="differential", 
            help="Fit type, defines POI and fit observable (recoil or mT)")
    parser.add_argument("--xsec", dest="xsec", action='store_true', 
            help="Write card for masked xsec normalization")
    parser.add_argument("--met", type=str, help="MET (DeepMETReso or RawPFMET)", default="RawPFMET")
    #parser.add_argument("--lumiScale", dest="lumiScale", help="Luminosity scale", type=float, default=1.0)
    return parser
    
def recoilSystNames(baseName, entries):
    systNames = []
    systNames.extend(["{baseName}_{i}{shift}".format(baseName=baseName, i=int(j/2), shift="Down" if j%2 else "Up") for j in range(entries)])
    return systNames

def main(args, xsec=False):
    logger = logging.setup_logger(__file__, args.verbose, args.noColorLogger)

    outfolder = f"CombineStudies/lowPU_{args.fitType}_{args.flavor}"
    if not os.path.isdir(outfolder):
        os.makedirs(outfolder)


    if not args.inputFile: args.inputFile = "lowPU_%s_%s.hdf5" % (args.flavor, args.met)

    datagroups = datagroupsLowPU(args.inputFile, flavor=args.flavor)
    unconstrainedProcs = [] # POI processes
    constrainedProcs = []   # constrained signal procs
    bkgProcs = ["Other", "Ztautau"] if args.fitType == "wmass" else ["Top", "EWK"] # background procs
    dataProc = "SingleMuon" if args.flavor == "mumu" else "SingleElectron"
    sigProc ="Zmumu" if args.flavor == "mumu" else "Zee"

    histName = "reco_mT"
    project = ["recoil_reco"]
    s = hist.tag.Slicer()
    if args.fitType == "differential":
        proc_base = dict(datagroups.groups[sigProc]) # signal process (Zmumu or Zee)
        for i in range(len(common.axis_recoil_gen_ptZ_lowpu)): # add gen bin processes to the dict
            proc_name = "Zmumu_genBin%d" % (i+1)
            proc_genbin = dict(proc_base)
            #proc_genbin['selectOp'] = lambda x, i=i: x[{"recoil_gen" : i}]
            proc_genbin['selectOpArgs'] = {"genBin": i}
            datagroups.addGroup(proc_name, proc_genbin)
            unconstrainedProcs.append(proc_name)
    elif args.fitType == "wmass": 
        constrainedProcs.append(sigProc)
        for proc in datagroups.groups.keys():
            datagroups.groups[proc]['selectOp'] = \
            lambda x: x[{ax : s[::hist.sum] for ax in ["recoil_gen", "mll",] if ax in x.axes.name}]
    elif args.fitType == "wlike":
        histName = "mT_corr_rec"
        project = ["mt"]
        constrainedProcs.append(sigProc) # need sum over gen bins
    elif args.fitType == "inclusive":
        unconstrainedProcs.append(sigProc)
        for proc in datagroups.groups.keys():
            datagroups.groups[proc]['selectOp'] = \
            lambda x, f=datagroups.groups[proc]['selectOp'] : f(x)[{ax : s[::hist.sum] for ax in ["recoil_gen", "recoil_reco",] if ax in x.axes.name}]

    
    suffix = ""
    if args.doStatOnly:
        suffix = "_stat"
    if xsec:
        suffix += "_xsec"
        bkgProcs = [] # for xsec norm card, remove all bkg procs but keep the data
        histName = "xnorm"
        project = ["count"]
        
        # fake data, as sum of all  Zmumu procs over recoil_gen
        proc_base = dict(datagroups.groups["Zmumu" if args.flavor == "mumu" else "Zee"])
        if args.fitType == "differential":
            proc_base['selectOp'] = lambda x, i=i: x[{"recoil_gen" : s[::hist.sum]}]
        dataProc = "fake_data"
        datagroups.addGroup(dataProc, proc_base)
    
    # hack: remove non-used procs/groups, as there can be more procs/groups defined than defined above
    # need to remove as cardTool takes all procs in the datagroups
    toDel = []
    for group in datagroups.groups: 
        if not group in constrainedProcs+unconstrainedProcs+bkgProcs+[dataProc]: toDel.append(group)
    datagroups.deleteGroup(toDel)    

    logger.debug(f"Going to use these groups: {datagroups.getNames()}")
    logger.debug(f"Datagroup keys: {datagroups.groups.keys()}")
    templateDir = f"{scriptdir}/Templates/LowPileupW"
    cardTool = CardTool.CardTool(f"{outfolder}/lowPU_{args.flavor}_{{chan}}_{args.met}_{args.fitType}{suffix}.txt")
    cardTool.setNominalTemplate(f"{templateDir}/main.txt")
    cardTool.setOutfile(os.path.abspath(f"{outfolder}/lowPU_{args.flavor}_{args.met}_{args.fitType}{suffix}.root"))
    cardTool.setProcesses(datagroups.getNames())
    cardTool.setDatagroups(datagroups)
    cardTool.setHistName(histName) 
    cardTool.setProjectionAxes(project) 
    cardTool.setNominalName(histName)
    #cardTool.setChannels([args.flavor])
    cardTool.setDataName(dataProc)
    cardTool.setProcsNoStatUnc(procs=[])
    cardTool.setSpacing(40)
    cardTool.setProcColumnsSpacing(20)
    cardTool.setWriteByCharge(True)
    cardTool.setUnconstrainedProcs(unconstrainedProcs)
    cardTool.setLumiScale(args.lumiScale)

    Zmumu_procs = cardTool.filteredProcesses(lambda x: "Zmumu" in x)
    Zmumu_procsIncTau = Zmumu_procs + ["Ztautau"]

    logger.debug(f"Making datacards with these processes: {cardTool.getProcesses()}")
    
    if args.fitType == "wlike" or args.fitType == "wmass":

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

    if args.doStatOnly:
        cardTool.addLnNSystematic("dummy", processes=cardTool.allMCProcesses(), size=1.001, group="dummy")
        cardTool.writeOutput(args=args)
        print("Using option --doStatOnly: the card was created with only mass weights and a dummy LnN syst on all processes")
        return

    pdfAction = {x : lambda h: h[{"recoil_gen" : s[::hist.sum]}] for x in Zmumu_procs if "gen" not in x},
    combine_helpers.add_pdf_uncertainty(cardTool, constrainedProcs+unconstrainedProcs, False, action=pdfAction)
    combine_helpers.add_scale_uncertainty(cardTool, args.minnlo_scale_unc, constrainedProcs+unconstrainedProcs, 
        to_fakes=False, use_hel_hist=True, resum=args.resumUnc)
    
    if not xsec:

        cardTool.addSystematic("prefireCorr",
            processes=cardTool.allMCProcesses(),
            mirror = False,
            group="CMS_prefire17",
            baseName="CMS_prefire17",
            systAxes = ["downUpVar"],
            labelsByAxis = ["downUpVar"],
        )
        
        combine_helpers.add_recoil_uncertainty(cardTool, Zmumu_procs, pu_type="lowPU")

        for lepEff in ["lepSF_HLT_DATA_stat", "lepSF_HLT_DATA_syst", "lepSF_HLT_MC_stat", "lepSF_HLT_MC_syst", "lepSF_ISO_stat", "lepSF_ISO_DATA_syst", "lepSF_ISO_MC_syst", "lepSF_IDIP_stat", "lepSF_IDIP_DATA_syst", "lepSF_IDIP_MC_syst"]:
            
            cardTool.addSystematic(lepEff,
                processes=cardTool.allMCProcesses(),
                mirror = True,
                group="CMS_lepton_eff",
                baseName=lepEff,
                systAxes = ["tensor_axis_0"],
                labelsByAxis = [""],
            )
        if args.fitType == "wmass":
            cardTool.addLnNSystematic("CMS_background", processes=["Other"], size=1.15, group="CMS_background")
        else:
            cardTool.addLnNSystematic("CMS_Top", processes=["Top"], size=1.06, group="CMS_background")
            cardTool.addLnNSystematic("CMS_VV", processes=["EWK"], size=1.16, group="CMS_background")

        cardTool.addLnNSystematic("CMS_lumi_lowPU", processes=cardTool.allMCProcesses(), size=1.02, group="CMS_lumi_lowPU")



    cardTool.writeOutput(args=args)

if __name__ == "__main__":
    parser = make_parser()
    args = parser.parse_args()

    main(args)
    if args.fitType == "differential":
        main(args, xsec=True)
        
    if args.fitType == "inclusive":
        main(args, xsec=True)    
        
        
