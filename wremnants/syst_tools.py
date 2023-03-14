import hist
import numpy as np
from utilities import boostHistHelpers as hh, common, logging
from wremnants import theory_tools
from wremnants.datasets.datagroups import datagroups2016
import collections.abc

logger = logging.child_logger(__name__)

def syst_transform_map(base_hist, hist_name):
    pdfInfo = theory_tools.pdfMapExtended 
    pdfNames = [pdfInfo[k]["name"] for k in pdfInfo.keys()]

    def pdfUnc(h, pdfName):
        key =  list(pdfInfo.keys())[list(pdfNames).index(pdfName)]
        unc = pdfInfo[key]["combine"]
        scale = pdfInfo[key]["scale"] if "scale" in pdfInfo[key] else 1.
        return theory_tools.hessianPdfUnc(h, uncType=unc, scale=scale)

    def uncHist(unc):
        return unc if base_hist == "nominal" else f"{base_hist}_{unc}"

    transforms = {}
    transforms.update({pdf+"Up" : {"action" : lambda h,p=pdf: pdfUnc(h, p)[0]} for pdf in pdfNames})
    transforms.update({pdf+"Down" : {"action" : lambda h,p=pdf: pdfUnc(h, p)[1]} for pdf in pdfNames})
    transforms.update({
        "massShift100MeVDown" : {"hist" : "massWeight", "action" : lambda h: h[{"tensor_axis_0" : 0}]},
        "massShift100MeVUp" : {"hist" : "massWeight", "action" : lambda h: h[{"tensor_axis_0" : 20}]},
    })

    s = hist.tag.Slicer()
    transforms.update({
        "QCDscale_muRmuFUp" : {
            "action" : lambda h: h[{"muRfact" : 2.j, "muFfact" : 2.j, "ptVgen" : s[::hist.sum]}]},
        "QCDscale_muRmuFDown" : {
            "action" : lambda h: h[{"muRfact" : 0.5j, "muFfact" : 0.5j, "ptVgen" : s[::hist.sum]}]},
        "QCDscale_muRUp" : {
            "action" : lambda h: h[{"muRfact" : 2.j, "muFfact" : 1.j, "ptVgen" : s[::hist.sum]}]},
        "QCDscale_muRDown" : {
            "action" : lambda h: h[{"muRfact" : 0.5j, "muFfact" : 1.j, "ptVgen" : s[::hist.sum]}]},
        "QCDscale_muFUp" : {
            "action" : lambda h: h[{"muRfact" : 1.j, "muFfact" : 2.j, "ptVgen" : s[::hist.sum]}]},
        "QCDscale_muFDown" : {
            "action" : lambda h: h[{"muRfact" : 1.j, "muFfact" : 0.5j, "ptVgen" : s[::hist.sum]}]},
        "QCDscale_cen" : {
            "action" : lambda h: h[{"muRfact" : 1.j, "muFfact" : 1.j, "ptVgen" : s[::hist.sum]}]},
    })

    def scetlibIdx(h, i):
        return h if not ("vars" in h.axes.name and h.axes["vars"].size > i) else h[{"vars" : i}]

    def projAx(hname):
        return hname.split("-")

    transforms.update({
        "resumFOScaleUp" : {
            "action" : lambda h: scetlibIdx(h, 2)},
        "resumFOScaleDown" : {
            "action" : lambda h: scetlibIdx(h, 1)},
        "resumLambdaDown" : {
            "action" : lambda h: scetlibIdx(h, 3)},
        "resumLambdaUp" : {
            "action" : lambda h: scetlibIdx(h, 4)},
        "resumTransitionUp" : {
            "action" : lambda h: hh.syst_min_or_max_env_hist(h, projAx(hist_name), "vars", 
                ["transition_points0.2_0.65_1.1", "transition_points0.4_0.55_0.7", 
                "transition_points0.2_0.45_0.7", "transition_points0.4_0.75_1.1", ],
                 no_flow=["ptVgen"], do_min=False)},
        "resumTransitionDown" : {
            "action" : lambda h: hh.syst_min_or_max_env_hist(h, projAx(hist_name), "vars", 
                ["transition_points0.2_0.65_1.1", "transition_points0.4_0.55_0.7", 
                "transition_points0.2_0.45_0.7", "transition_points0.4_0.75_1.1", ],
                 no_flow=["ptVgen"], do_min=True)},
       "resumNPUp" : {
            "action" : lambda h: hh.syst_min_or_max_env_hist(h, projAx(hist_name), "vars", 
                ['c_nu-0.15-omega_nu0.43', 'c_nu0.05', 'c_nu0.5-omega_nu0.15', 'c_nu-0.5-omega_nu0.37'],
                 no_flow=["ptVgen"], do_min=False)},
        "resumNPDown" : {
            "action" : lambda h: hh.syst_min_or_max_env_hist(h, projAx(hist_name), "vars", 
                ['c_nu-0.15-omega_nu0.43', 'c_nu0.05', 'c_nu0.5-omega_nu0.15', 'c_nu-0.5-omega_nu0.37'],
                 no_flow=["ptVgen"], do_min=True)},
        "resumScaleMax" : {
            "action" : lambda h: hh.syst_min_or_max_env_hist(h, projAx(hist_name), "vars", range(9,44), no_flow=["ptVgen"], do_min=False)},
        "resumScaleMin" : {
            "action" : lambda h: hh.syst_min_or_max_env_hist(h, projAx(hist_name), "vars", range(9,44), no_flow=["ptVgen"], do_min=True)},
    })
    for k,v in transforms.items():
        if any([x in k for x in ["QCDscale", "resum", "pdf"]]):
            v["procs"] = common.vprocs 

    return transforms

def scale_helicity_hist_to_variations(scale_hist, sum_axes=[], rebinPtV=None):
    
    s = hist.tag.Slicer()
    # select nominal QCD scales, but keep the sliced axis at size 1 for broadcasting
    nom_scale_hist = scale_hist[{"muRfact" : s[1.j:1.j+1], "muFfact" : s[1.j:1.j+1]}]
    axisNames = scale_hist.axes.name
    # select nominal QCD scales and project down to nominal axes
    genAxes = ["ptVgen", "chargeVgen", "helicity"]
    nom_hist = nom_scale_hist[{"muRfact" : s[1.j], "muFfact" : s[1.j] }]
    for genAxis in genAxes:
        if genAxis in axisNames:
            nom_hist = nom_hist[{genAxis : s[::hist.sum]}]
    
    hasHelicityAxis = "helicity" in axisNames
    hasPtAxis = "ptVgen" in axisNames

    if rebinPtV and hasPtAxis:
        # Treat single bin array as a float
        array_rebin = isinstance(rebinPtV, collections.abc.Sequence)
        if array_rebin and len(rebinPtV) == 1:
            rebinPtV = rebinPtV[0]
            array_rebin = False

        if array_rebin:
            scale_hist = hh.rebinHist(scale_hist, "ptVgen", rebinPtV)
            nom_scale_hist = hh.rebinHist(nom_scale_hist, "ptVgen", rebinPtV)
        else:
            scale_hist = scale_hist[{"ptVgen" : s[::hist.rebin(rebinPtV)]}]
            nom_scale_hist = nom_scale_hist[{"ptVgen" : s[::hist.rebin(rebinPtV)]}]

    for axis in sum_axes:
        if axis in axisNames:
            scale_hist = scale_hist[{axis : s[::hist.sum]}]
            nom_scale_hist = nom_scale_hist[{axis : s[::hist.sum]}]
        else:
            logger.warning(f"In scale_helicity_hist_to_variations: axis '{axis}' not found in histogram.")
        
    # difference between a given scale and the nominal, plus the sum
    # this emulates the "weight if idx else nominal" logic and corresponds to the decorrelated
    # variations
    if scale_hist.name is None:
        out_name = "scale_helicity_variations" if hasHelicityAxis else "scale_vpt_variations" if hasPtAxis else "scale_vcharge_variations"
    else:
        out_name = scale_hist.name + "_variations"

    expd = scale_hist.ndim - nom_hist.ndim
    expandnom = np.expand_dims(nom_hist.view(flow=True), [-expd+i for i in range(expd)])
    systhist = scale_hist.view(flow=True) - nom_scale_hist.view(flow=True) + expandnom

    scale_variation_hist = hist.Hist(*scale_hist.axes, storage = scale_hist._storage_type(), 
                                     name = out_name, data = systhist)

    return scale_variation_hist 

def uncertainty_hist_from_envelope(h, proj_ax, entries):
    hdown = hh.syst_min_or_max_env_hist(h, proj_ax, "vars", entries, no_flow=["ptVgen"], do_min=True)
    hup = hh.syst_min_or_max_env_hist(h, proj_ax, "vars", entries, no_flow=["ptVgen"], do_min=False)
    hnew = hist.Hist(*h.axes[:-1], common.down_up_axis, storage=h._storage_type())
    hnew[...,0] = hdown.view(flow=True)
    hnew[...,1] = hup.view(flow=True)
    return hnew

def define_mass_weights(df, proc):
    nweights = 23 if proc in common.zprocs_all else 21
    # from -100 to 100 MeV with 10 MeV increment
    df = df.Define("massWeight_tensor", f"wrem::vec_to_tensor_t<double, {nweights}>(MEParamWeight)")
    df = df.Define("massWeight_tensor_wnom", "auto res = massWeight_tensor; res = nominal_weight*res; return res;")

    return df

def add_massweights_hist(results, df, axes, cols, base_name="nominal", proc=""):
    name = datagroups2016.histName(base_name, syst="massWeight")
    massWeight = df.HistoBoost(name, axes, [*cols, "massWeight_tensor_wnom"], 
                    tensor_axes=[hist.axis.StrCategory(massWeightNames(proc=proc), name="massShift")])
    results.append(massWeight)

def massWeightNames(matches=None, proc=""):
    central=10
    nweights=21
    names = [f"massShift{int(abs(central-i)*10)}MeV{'' if i == central else ('Down' if i < central else 'Up')}" for i in range(nweights)]
    if proc and proc in common.zprocs_all:
        # This is the PDG uncertainty (turned off for now since it doesn't seem to have been read into the nano)
        names.extend(["massShift2p1MeVDown", "massShift2p1MeVUp"])

    # If name is "" it won't be stored
    return [x if not matches or any(y in x for y in matches) else "" for x in names]

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

        name = datagroups2016.histName(base_name, syst=pdfName)
        names = [f"pdf{pdfName}{i}" for i in range(pdfInfo["entries"])] if pdfInfo["combine"] == "symHessian" else \
                theory_tools.pdfNamesAsymHessian(pdfInfo["entries"], pdfName)
        pdf_ax = hist.axis.StrCategory(names, name="pdfVar")
        pdfHist = df.HistoBoost(name, axes, [*cols, tensorName], tensor_axes=[pdf_ax])

        if pdfInfo["alphasRange"] == "001":
            name = datagroups2016.histName(base_name, syst=f"{pdfName}alphaS001")
            as_ax = hist.axis.StrCategory(["as0117", "as0119"], name="alphasVar")
        else:
            name = datagroups2016.histName(base_name, syst=f"{pdfName}alphaS002")
            as_ax = hist.axis.StrCategory(["as0116", "as0120"], name="alphasVar")
        alphaSHist = df.HistoBoost(name, axes, [*cols, tensorASName], tensor_axes=[as_ax])
        results.extend([pdfHist, alphaSHist])
    return df

def add_qcdScale_hist(results, df, axes, cols, base_name="nominal"):
    name = datagroups2016.histName(base_name, syst="qcdScale")
    scaleHist = df.HistoBoost(name, axes, [*cols, "scaleWeights_tensor_wnom"], tensor_axes=theory_tools.scale_tensor_axes)
    results.append(scaleHist)

def add_qcdScaleByHelicityUnc_hist(results, df, helper, axes, cols, base_name="nominal"):
    name = datagroups2016.histName(base_name, syst="qcdScaleByHelicity")
    df = df.Define("helicityWeight_tensor", helper, ["massVgen", "absYVgen", "ptVgen", "chargeVgen", "csSineCosThetaPhi", "scaleWeights_tensor", "nominal_weight"])
    qcdScaleByHelicityUnc = df.HistoBoost(name, axes, [*cols,"helicityWeight_tensor"], tensor_axes=helper.tensor_axes)
    results.append(qcdScaleByHelicityUnc)

def add_QCDbkg_jetPt_hist(results, df, nominal_axes, nominal_cols, base_name="nominal", jet_pt=30):
    # branching the rdataframe to add special filter, no need to return dQCDbkGVar
    name = datagroups2016.histName(base_name, syst=f"qcdJetPt{str(jet_pt)}")
    dQCDbkGVar = df.Define(f"goodCleanJetsPt{jet_pt}", f"goodCleanJetsNoPt && Jet_pt > {jet_pt}")
    dQCDbkGVar = dQCDbkGVar.Filter(f"passMT || Sum(goodCleanJetsPt{jet_pt})>=1")
    qcdJetPt = dQCDbkGVar.HistoBoost(name, nominal_axes, [*nominal_cols, "nominal_weight"])
    results.append(qcdJetPt)
                                        
def add_muon_efficiency_unc_hists(results, df, helper_stat, helper_syst, axes, cols, base_name="nominal", is_w_like=False):

    if is_w_like:
        muon_columns_stat = ["trigMuons_pt0", "trigMuons_eta0", "trigMuons_charge0", "nonTrigMuons_pt0", "nonTrigMuons_eta0", "nonTrigMuons_charge0"]
        muon_columns_syst = ["trigMuons_pt0", "trigMuons_eta0", "trigMuons_SApt0", "trigMuons_SAeta0", "trigMuons_charge0",
            "nonTrigMuons_pt0", "nonTrigMuons_eta0", "nonTrigMuons_SApt0", "nonTrigMuons_SAeta0", "nonTrigMuons_charge0"]
    else:
        muon_columns_stat = ["goodMuons_pt0", "goodMuons_eta0", "goodMuons_charge0"]
        muon_columns_syst = ["goodMuons_pt0", "goodMuons_eta0", "goodMuons_SApt0", "goodMuons_SAeta0", "goodMuons_charge0", "passIso"]

    for key,helper in helper_stat.items():
        if "iso" in key and not is_w_like:
            df = df.Define(f"effStatTnP_{key}_tensor", helper, [*muon_columns_stat, "passIso", "nominal_weight"])        
        else:
            df = df.Define(f"effStatTnP_{key}_tensor", helper, [*muon_columns_stat, "nominal_weight"])
        name = datagroups2016.histName(base_name, syst=f"effStatTnP_{key}")
        effStatTnP = df.HistoBoost(name, axes, [*cols, f"effStatTnP_{key}_tensor"], tensor_axes = helper.tensor_axes)
        results.append(effStatTnP)
    
    df = df.Define("effSystTnP_weight", helper_syst, [*muon_columns_syst, "nominal_weight"])
    name = datagroups2016.histName(base_name, syst=f"effSystTnP")
    effSystTnP = df.HistoBoost(name, axes, [*cols, "effSystTnP_weight"], tensor_axes = helper_syst.tensor_axes)
    results.append(effSystTnP)

    return df

def add_L1Prefire_unc_hists(results, df, helper_stat, helper_syst, axes, cols, base_name="nominal"):
    df = df.Define("muonL1PrefireStat_tensor", helper_stat, ["Muon_correctedEta", "Muon_correctedPt", "Muon_correctedPhi", "Muon_correctedCharge", "Muon_looseId", "nominal_weight"])
    name = datagroups2016.histName(base_name, syst=f"muonL1PrefireStat")
    muonL1PrefireStat = df.HistoBoost(name, axes, [*cols, "muonL1PrefireStat_tensor"], tensor_axes = helper_stat.tensor_axes)
    results.append(muonL1PrefireStat)

    df = df.Define("muonL1PrefireSyst_tensor", helper_syst, ["Muon_correctedEta", "Muon_correctedPt", "Muon_correctedPhi", "Muon_correctedCharge", "Muon_looseId", "nominal_weight"])
    name = datagroups2016.histName(base_name, syst=f"muonL1PrefireSyst")
    muonL1PrefireSyst = df.HistoBoost(name, axes, [*cols, "muonL1PrefireSyst_tensor"], tensor_axes = [common.down_up_axis])
    results.append(muonL1PrefireSyst)

    df = df.Define("ecalL1Prefire_tensor", f"wrem::twoPointScaling(nominal_weight/L1PreFiringWeight_ECAL_Nom, L1PreFiringWeight_ECAL_Dn, L1PreFiringWeight_ECAL_Up)")
    name = datagroups2016.histName(base_name, syst=f"ecalL1Prefire")
    ecalL1Prefire = df.HistoBoost(name, axes, [*cols, "ecalL1Prefire_tensor"], tensor_axes = [common.down_up_axis])
    results.append(ecalL1Prefire)

    return df

def add_muonscale_hist(results, df, netabins, mag, isW, axes, cols, base_name="nominal", muon_eta="goodMuons_eta0"):
    nweights = 21 if isW else 23

    df = df.Define(f"muonScaleDummy{netabins}Bins{muon_eta}", f"wrem::dummyScaleFromMassWeights<{netabins}, {nweights}>(nominal_weight, massWeight_tensor, {muon_eta}, {mag}, {str(isW).lower()})")

    scale_etabins_axis = hist.axis.Regular(netabins, -2.4, 2.4, name="scaleEtaSlice", underflow=False, overflow=False)
    name = datagroups2016.histName(base_name, syst=f"muonScaleSyst")

    dummyMuonScaleSyst = df.HistoBoost(name, axes, [*cols, f"muonScaleDummy{netabins}Bins{muon_eta}"], tensor_axes=[common.down_up_axis, scale_etabins_axis])
    results.append(dummyMuonScaleSyst)

    return df


def add_muonscale_smeared_hist(results, df, netabins, mag, isW, axes, cols, base_name="nominal", muon_eta="goodMuons_eta0"):
    # add_muonscale_hist has to be called first such that "muonScaleDummy{netabins}Bins{muon_eta}" is defined
    nweights = 21 if isW else 23

    scale_etabins_axis = hist.axis.Regular(netabins, -2.4, 2.4, name="scaleEtaSlice", underflow=False, overflow=False)
    name = datagroups2016.histName(base_name, syst=f"muonScaleSyst_gen_smear")

    dummyMuonScaleSyst_gen_smear = df.HistoBoost(name, axes, [*cols, f"muonScaleDummy{netabins}Bins{muon_eta}"], tensor_axes=[common.down_up_axis, scale_etabins_axis])
    results.append(dummyMuonScaleSyst_gen_smear)

    return df

def scetlib_scale_unc_hist(h, obs, syst_ax="vars"):
    hnew = hist.Hist(*h.axes[:-1], hist.axis.StrCategory(["central"]+scetlib_scale_vars(),
    					name=syst_ax), storage=h._storage_type())
    
    hnew[...,"central"] = h[...,"central"].view(flow=True)
    hnew[...,"resumFOScaleUp"] = h[...,"kappaFO2."].view(flow=True)
    hnew[...,"resumFOScaleDown"] = h[...,"kappaFO0.5"].view(flow=True)
    hnew[...,"resumLambdaUp"] = h[...,"lambda0.8"].view(flow=True)
    hnew[...,"resumLambdaDown"] = h[...,"lambda1.5"].view(flow=True)
    
    transition_names = [x for x in h.axes[syst_ax] if "transition" in x]    
    hnew[...,"resumTransitionUp"] = hh.syst_min_or_max_env_hist(h, obs, syst_ax, 
                                    h.axes[syst_ax].index(transition_names), do_min=False).view(flow=True)
    hnew[...,"resumTransitionDown"] = hh.syst_min_or_max_env_hist(h, obs, syst_ax, 
                                    h.axes[syst_ax].index(transition_names), do_min=True).view(flow=True)
    
    resum_names = [x for x in h.axes[syst_ax] if not any(i in x for i in ["lambda", "kappa", "transition"])]
    hnew[...,"resumScaleUp"] = hh.syst_min_or_max_env_hist(h, obs, syst_ax, 
                                    h.axes[syst_ax].index(resum_names), do_min=False).view(flow=True)
    hnew[...,"resumScaleDown"] = hh.syst_min_or_max_env_hist(h, obs, syst_ax, 
                                    h.axes[syst_ax].index(resum_names), do_min=True).view(flow=True)
    return hnew
