import hist
import numpy as np
from utilities import boostHistHelpers as hh, common, logging
from wremnants import theory_tools, syst_tools
from wremnants.datasets.datagroups import Datagroups
from wremnants.helicity_utils import *
import re
import collections.abc

logger = logging.child_logger(__name__)


def add_massweights_hist(results, df, axes, cols, base_name="nominal", proc=""):
    name = Datagroups.histName(base_name, syst="massWeight")
    mass_axis = hist.axis.StrCategory(syst_tools.massWeightNames(proc=proc), name="massShift")
    massweightHelicity, massWeight_axes = make_massweight_helper_helicity(mass_axis)
    df=df.Define("massWeight_tensor_wnom_helicity", massweightHelicity, ['massWeight_tensor_wnom', 'helWeight_tensor'])
    massWeight = df.HistoBoost(name, axes, [*cols, "massWeight_tensor_wnom_helicity"], 
                    tensor_axes=massWeight_axes, 
                    storage=hist.storage.Double())
    results.append(massWeight)

def add_pdf_hists(results, df, dataset, axes, cols, pdfs, base_name="nominal"):
    for pdf in pdfs:
        try:
            pdfInfo = theory_tools.pdf_info_map(dataset, pdf)
        except ValueError as e:
            logger.info(e)
            continue

        pdfName = pdfInfo["name"]
        tensorName = f"{pdfName}Weights_tensor"
        tensorASName = f"{pdfName}ASWeights_tensor"
        npdf=pdfInfo["entries"]
        name = Datagroups.histName(base_name, syst=pdfName)
        names = getattr(theory_tools, f"pdfNames{'Sym' if pdfInfo['combine'] == 'symHessian' else 'Asym'}Hessian")(pdfInfo["entries"], pdfName)
        pdf_ax = hist.axis.StrCategory(names, name="pdfVar")
        if tensorName not in df.GetColumnNames():
            logger.warning(f"PDF {pdf} was not found for sample {dataset}. Skipping uncertainty hist!")
            continue
        pdfHeltensor, pdfHeltensor_axes =  make_pdfweight_helper_helicity(npdf, pdf_ax)
        df = df.Define(f'{tensorName}_helicity', pdfHeltensor, [tensorName, "helWeight_tensor"])
        pdfHist = df.HistoBoost(name, axes, [*cols, f'{tensorName}_helicity'], tensor_axes=pdfHeltensor_axes, storage=hist.storage.Double())

        if pdfInfo["alphasRange"] == "001":
            name = Datagroups.histName(base_name, syst=f"{pdfName}alphaS001")
            as_ax = hist.axis.StrCategory(["as0117", "as0119"], name="alphasVar")
        else:
            name = Datagroups.histName(base_name, syst=f"{pdfName}alphaS002")
            as_ax = hist.axis.StrCategory(["as0116", "as0120"], name="alphasVar")
        alphaSHeltensor, alphaSHeltensor_axes =  make_pdfweight_helper_helicity(2, as_ax)
        df = df.Define(f'{tensorASName}_helicity', alphaSHeltensor, [tensorASName, "helWeight_tensor"])
        alphaSHist = df.HistoBoost(name, axes, [*cols, f'{tensorASName}_helicity'], tensor_axes=alphaSHeltensor_axes, storage=hist.storage.Double())
        results.extend([pdfHist, alphaSHist])
    return df

def add_qcdScale_hist(results, df, axes, cols, base_name="nominal"):
    name = Datagroups.histName(base_name, syst="qcdScale")
    qcdbyHelicity, qcdbyHelicity_axes = make_qcdscale_helper_helicity(theory_tools.scale_tensor_axes)
    df = df.Define('scaleWeights_tensor_wnom_helicity', qcdbyHelicity, ['scaleWeights_tensor_wnom', 'helWeight_tensor'])
    scaleHist = df.HistoBoost(name, axes, [*cols, "scaleWeights_tensor_wnom_helicity"], tensor_axes=qcdbyHelicity_axes, storage=hist.storage.Double())
    results.append(scaleHist)

def add_muon_efficiency_unc_hists(results, df, helper_stat, helper_syst, axes, cols, base_name="nominal", is_w_like=False):
    if is_w_like:
        muon_columns_stat = ["trigMuons_pt0", "trigMuons_eta0", "trigMuons_charge0", "nonTrigMuons_pt0", "nonTrigMuons_eta0", "nonTrigMuons_charge0"]
        muon_columns_syst = ["trigMuons_pt0", "trigMuons_eta0", "trigMuons_SApt0", "trigMuons_SAeta0", "trigMuons_charge0",
            "nonTrigMuons_pt0", "nonTrigMuons_eta0", "nonTrigMuons_SApt0", "nonTrigMuons_SAeta0", "nonTrigMuons_charge0"]
    else:
        # FIXME: this should read standalone variables when effStat is for tracking
        muon_columns_stat = ["goodMuons_pt0", "goodMuons_eta0", "goodMuons_charge0"]
        muon_columns_syst = ["goodMuons_pt0", "goodMuons_eta0", "goodMuons_SApt0", "goodMuons_SAeta0", "goodMuons_charge0", "passIso"]

    for key,helper in helper_stat.items():
        if "iso" in key and not is_w_like:
            df = df.Define(f"effStatTnP_{key}_tensor", helper, [*muon_columns_stat, "passIso", "nominal_weight"])        
        else:
            df = df.Define(f"effStatTnP_{key}_tensor", helper, [*muon_columns_stat, "nominal_weight"])
        helper_helicity, helper_helicity_axes = make_muon_eff_stat_helpers_helicity(helper)
        df = df.Define(f"effStatTnP_{key}_ByHelicity_tensor", helper_helicity, [f"effStatTnP_{key}_tensor", "helWeight_tensor"])
        name = Datagroups.histName(base_name, syst=f"effStatTnP_{key}")
        effStatTnP = df.HistoBoost(name, axes, [*cols, f"effStatTnP_{key}_ByHelicity_tensor"], tensor_axes = helper_helicity_axes, storage=hist.storage.Double())
        results.append(effStatTnP)
    
    df = df.Define("effSystTnP_weight", helper_syst, [*muon_columns_syst, "nominal_weight"])

    helper_syst_helicity, helper_syst_helicity_axes = make_muon_eff_syst_helper_helicity(helper_syst)
    df = df.Define("effSystTnP_weight_ByHelicity_tensor", helper_syst_helicity, ["effSystTnP_weight", "helWeight_tensor"])
    name = Datagroups.histName(base_name, syst=f"effSystTnP")
    effSystTnP = df.HistoBoost(name, axes, [*cols, "effSystTnP_weight_ByHelicity_tensor"], tensor_axes = helper_syst_helicity_axes)
    results.append(effSystTnP)
    
    return df

def add_L1Prefire_unc_hists(results, df, helper_stat, helper_syst, axes, cols, base_name="nominal"):
    df = df.Define("muonL1PrefireStat_tensor", helper_stat, ["Muon_correctedEta", "Muon_correctedPt", "Muon_correctedPhi", "Muon_correctedCharge", "Muon_looseId", "nominal_weight"])
    prefirebyhelicity_stat, prefire_axes_stat = make_muon_prefiring_helper_stat_byHelicity(helper_stat)
    df = df.Define("muonL1PrefireStatByHelicity_tensor", prefirebyhelicity_stat, ["muonL1PrefireStat_tensor", "helWeight_tensor"])
    name = Datagroups.histName(base_name, syst=f"muonL1PrefireStat")
    muonL1PrefireStat = df.HistoBoost(name, axes, [*cols, "muonL1PrefireStatByHelicity_tensor"], tensor_axes = prefire_axes_stat, storage=hist.storage.Double())
    results.append(muonL1PrefireStat)

    df = df.Define("muonL1PrefireSyst_tensor", helper_syst, ["Muon_correctedEta", "Muon_correctedPt", "Muon_correctedPhi", "Muon_correctedCharge", "Muon_looseId", "nominal_weight"])
    prefirebyhelicity_syst, prefire_axes_syst = make_muon_prefiring_helper_syst_byHelicity()
    df = df.Define("muonL1PrefireSystByHelicity_tensor", prefirebyhelicity_syst, ["muonL1PrefireSyst_tensor", "helWeight_tensor"])
    name = Datagroups.histName(base_name, syst=f"muonL1PrefireSyst")
    muonL1PrefireSyst = df.HistoBoost(name, axes, [*cols, "muonL1PrefireSystByHelicity_tensor"], tensor_axes = prefire_axes_syst, storage=hist.storage.Double())
    results.append(muonL1PrefireSyst)

    df = df.Define("ecalL1Prefire_tensor", f"wrem::twoPointScaling(nominal_weight/L1PreFiringWeight_ECAL_Nom, L1PreFiringWeight_ECAL_Dn, L1PreFiringWeight_ECAL_Up)")
    #can reuse the same helper since it's the tensor multiplication of same types
    df = df.Define("ecalL1PrefireByHelicity_tensor", prefirebyhelicity_syst, ["ecalL1Prefire_tensor", "helWeight_tensor"])
    name = Datagroups.histName(base_name, syst=f"ecalL1Prefire")
    ecalL1Prefire = df.HistoBoost(name, axes, [*cols, "ecalL1PrefireByHelicity_tensor"], tensor_axes = prefire_axes_syst, storage=hist.storage.Double())
    results.append(ecalL1Prefire)

    return df
'''
def add_muonscale_hist(results, df, netabins, mag, isW, axes, cols, base_name="nominal", muon_eta="goodMuons_eta0"):
    nweights = 21 if isW else 23

    df = df.Define(f"muonScaleDummy{netabins}Bins{muon_eta}", f"wrem::dummyScaleFromMassWeights<{netabins}, {nweights}>(nominal_weight, massWeight_tensor, {muon_eta}, {mag}, {str(isW).lower()})")

    scale_etabins_axis = hist.axis.Regular(netabins, -2.4, 2.4, name="scaleEtaSlice", underflow=False, overflow=False)
    name = Datagroups.histName(base_name, syst=f"muonScaleSyst")

    dummyMuonScaleSyst = df.HistoBoost(name, axes, [*cols, f"muonScaleDummy{netabins}Bins{muon_eta}"], tensor_axes=[common.down_up_axis, scale_etabins_axis], storage=hist.storage.Double())
    results.append(dummyMuonScaleSyst)

    return df


def add_muonscale_smeared_hist(results, df, netabins, mag, isW, axes, cols, base_name="nominal", muon_eta="goodMuons_eta0"):
    # add_muonscale_hist has to be called first such that "muonScaleDummy{netabins}Bins{muon_eta}" is defined
    nweights = 21 if isW else 23

    scale_etabins_axis = hist.axis.Regular(netabins, -2.4, 2.4, name="scaleEtaSlice", underflow=False, overflow=False)
    name = Datagroups.histName(base_name, syst=f"muonScaleSyst_gen_smear")

    dummyMuonScaleSyst_gen_smear = df.HistoBoost(name, axes, [*cols, f"muonScaleDummy{netabins}Bins{muon_eta}"], tensor_axes=[common.down_up_axis, scale_etabins_axis], storage=hist.storage.Double())
    results.append(dummyMuonScaleSyst_gen_smear)

    return df
'''
def add_theory_hists(results, df, args, dataset_name, corr_helpers, qcdScaleByHelicity_helper, axes, cols, base_name="nominal", for_wmass=True):
    axis_chargeVgen = qcdScaleByHelicity_helper.hist.axes["chargeVgen"]
    axis_ptVgen = hist.axis.Variable(
        common.ptV_10quantiles_binning, 
        name = "ptVgen", underflow=False
    )
    #difference from the other, since qt should already be in the cols
    scale_axes = [*axes, axis_chargeVgen]
    scale_cols = [*cols, "chargeVgen"]

    df = theory_tools.define_scale_tensor(df)
    df = syst_tools.define_mass_weights(df, dataset_name)

    add_pdf_hists(results, df, dataset_name, axes, cols, args.pdfs, base_name=base_name)
    add_qcdScale_hist(results, df, scale_axes, scale_cols, base_name=base_name)

    #isZ = dataset_name in common.zprocs

    #if args.theoryCorr and dataset_name in corr_helpers:
    #    results.extend(theory_tools.make_theory_corr_hists(df, base_name, axes, cols, 
    #        corr_helpers[dataset_name], args.theoryCorr, modify_central_weight=not args.theoryCorrAltOnly, isW = not isZ)
    #    )

    #if for_wmass or isZ:
        # there is no W backgrounds for the Wlike, make QCD scale histograms only for Z
        # should probably remove the charge here, because the Z only has a single charge and the pt distribution does not depend on which charged lepton is selected


        # TODO: Should have consistent order here with the scetlib correction function
    add_massweights_hist(results, df, axes, cols, proc=dataset_name, base_name=base_name)

    return df

