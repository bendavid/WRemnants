from utilities import boostHistHelpers as hh, common, logging, input_tools
from wremnants import syst_tools,theory_tools,recoil_tools
import numpy as np
import re

logger = logging.child_logger(__name__)

def add_modeling_uncertainty(card_tool, minnlo_scale, signal_samples, background_samples, to_fakes, resumType, wmass, scaleTNP=1, rebin_pt=None):
    scale_name = "W" if wmass else "Z"
    resum_samples = signal_samples+background_samples if wmass and resumType == "tnp" else signal_samples
    scale_label = scale_name if resumType != "tnp" else ""

    do_resum = resumType and resumType not in ["none", "npOnly"]
    resum_np = resumType and resumType != ["none"]

    if resumType in ["tnp", "npOnly", "scale"]:
        syst_ax = "vars"
        corr_hist_name = theory_unc_hist(card_tool)
        if "scetlib_dyturbo" not in corr_hist_name:
            raise ValueError("Resummation uncertainties only defined for SCETlib+DYturbo")

        corr_hist = card_tool.getHistsForProcAndSyst(signal_samples[0], corr_hist_name)


    if do_resum:
        add_resum_uncertainty(card_tool, resum_samples, to_fakes, 
                                            uncType=resumType, scale=scaleTNP, name_append=scale_label)
        # Don't correlate transition uncertainty
        add_resum_transition_uncertainty(card_tool, signal_samples, to_fakes, name_append=scale_name)

    add_minnlo_scale_uncertainty(card_tool, minnlo_scale, signal_samples, to_fakes, 
                                          resum=do_resum, name_append=scale_name)
    # for Z background in W mass case (W background for Wlike is essentially 0, useless to apply QCD scales there)
    # For now the background (Z) is always MiNNLO scale uncertainty
    if wmass and background_samples:
        add_minnlo_scale_uncertainty(card_tool, minnlo_scale, background_samples, to_fakes, 
                                            resum=do_resum, name_append="Z")
        add_resum_transition_uncertainty(card_tool, background_samples, to_fakes, name_append="Z")

        if resumType != "tnp":
            add_resum_uncertainty(card_tool, background_samples, to_fakes, 
                                            uncType=resumType, scale=scaleTNP, name_append="Z")

    if resumType != "none":
        add_common_np_uncertainties(card_tool, signal_samples+background_samples, to_fakes)

        omega = False
        npfunc = add_OmegaNP_uncertainties if omega else add_LambdaNP_uncertainties
        npfunc(card_tool, signal_samples, to_fakes, name_append=scale_name)

        if wmass and background_samples:
            npfunc(card_tool, background_samples, to_fakes, name_append="Z")

def match_str_axis_entries(str_axis, match_re):
    return [x for x in str_axis if any(re.match(r, x) for r in match_re)]

def add_resum_uncertainty(card_tool, samples, to_fakes, uncType, name_append="", scale=1):
    if input_tools.args_from_metadata(card_tool, "theoryCorrAltOnly"):
        logger.error("The theory correction was only applied as an alternate hist. Using its syst isn't well defined!")


    syst_ax = "vars"
    np_nuisances = ["^c_nu-*\d+", "^omega_nu-*\d+", "^Omega-*\d+"]
    both_exclude = ['^kappaFO.*','^recoil_scheme.*',"^transition_points.*",]+np_nuisances
    tnp_nuisance_names = ['gamma_cusp', 'gamma_mu_q', 'gamma_nu', 'h_qqV', 's', 'b_qqV', 'b_qqbarV', 'b_qqS', 'b_qqDS', 'b_qg', 'b_qg', ]

def add_resum_scale_uncertainty(card_tool, theory_hist_name, samples, to_fakes, name_append=""):
    obs = card_tool.project[:]
    if not obs:
        raise ValueError("Failed to find the observable names for the resummation uncertainties")
    
    theory_hist = card_tool.getHistsForProcAndSyst(samples[0], theory_hist_name)
    resumscale_nuisances = match_str_axis_entries(h.axes[syst_ax], ["^nuB.*", "nuS.*", "^muB.*", "^muS.*",])

    expanded_samples = card_tool.datagroups.getProcNames(samples)
    syst_ax = "vars"

    card_tool.addSystematic(name=theory_hist,
        processes=samples,
        group="resumScale",
        passToFakes=to_fakes,
        skipEntries=[{syst_ax : x} for x in both_exclude+tnp_nuisances],
        systAxes=["downUpVar"], # Is added by the actionMap
        actionMap={s : lambda h: hh.syst_min_and_max_env_hist(h, obs, "vars", resumscale_nuisances) for s in expanded_samples},
        outNames=[f"scetlibResumScale{name_append}Up", f"scetlibResumScale{name_append}Down"],
        rename=f"resumScale{name_append}",
        systNamePrepend=f"resumScale{name_append}_",
    )
    #TODO: check if this is actually the proper treatment of these uncertainties
    card_tool.addSystematic(name=theory_hist,
        processes=samples,
        group="resumScale",
        passToFakes=to_fakes,
        systAxes=["vars"],
        actionMap={s : lambda h: h[{"vars" : ["kappaFO0.5-kappaf2.", "kappaFO2.-kappaf0.5", "mufdown", "mufup",]}] for s in expanded_samples},
        outNames=[f"scetlib_kappa{name_append}Up", f"scetlib_kappa{name_append}Down", f"scetlib_muF{name_append}Up", f"scetlib_muF{name_append}Down"],
        rename=f"resumFOScale{name_append}",
        systNamePrepend=f"resumScale{name_append}_",
    )

def add_resum_transition_uncertainty(card_tool, samples, to_fakes, name_append=""):
    obs = card_tool.project[:]
    theory_hist = theory_unc_hist(card_tool)
    expanded_samples = card_tool.datagroups.getProcNames(samples)

    card_tool.addSystematic(name=theory_hist,
        processes=samples,
        group="resumTransition",
        systAxes=["downUpVar"],
        passToFakes=to_fakes,
        # NOTE: I don't actually remember why this used no_flow=ptVgen previously, I don't think there's any harm in not using it...
        actionMap={s : lambda h: hh.syst_min_and_max_env_hist(h, obs, "vars", 
            [x for x in h.axes["vars"] if "transition_point" in x]) for s in expanded_samples},
        outNames=[f"resumTransition{name_append}Up", f"resumTransition{name_append}Down"],
        rename=f"scetlibResumTransition{name_append}",
    )

def add_OmegaNP_uncertainties(card_tool, samples, to_fakes, name_append, mode="delta"):
    if mode not in ["binned", "unbinned", "delta"]:
        raise ValueError(f"mode '{mode}' is not a valid choice for OmegaNP parameters")

    theory_hist = theory_unc_hist(card_tool, "Omega" if mode == "binned" else "")

    omega_vals = ["Omega0.", "Omega0.8"]
    action = lambda h: h[{"vars" : omega_vals}]
    if mode == "binned":
        action=lambda h: syst_tools.hist_to_variations(h[{"vars" : ["central" if "central" in h.axes["vars"] else "pdf0", *omega_vals]}]),

    card_tool.addSystematic(name=hist,
        processes=samples,
        group="resumNonpert",
        systAxes=["vars",] if mode != "binned" else ["absYVgenNP", "chargeVgenNP", "vars"],
        passToFakes=to_fakes,
        action=action,
        # Careful here... it's nasty, but the order of the replace matters
        systNameReplace=[(omega_vals[1], "OmegaUp"), (omega_vals[0], "OmegaDown"), ],
        rename=f"scetlibOmegaNP{name_append}",
    )

    if delta:
        domega_name = f"scetlibDeltaOmegaNP{name_append}"
        card_tool.addSystematic(name=hist,
            processes=samples,
            group="resumNonpert",
            systAxes=["vars",],
            passToFakes=to_fakes,
            action=lambda h: h[{"vars" : ["Delta_Omega-0.02", "Delta_Omega0.02"]}],
            outNames=[domega_name+"Down", domega_name+"Up"],
            rename=domega_name,
        )

def add_LambdaNP_uncertainties(card_tool, samples, to_fakes, name_append):
    theory_hist = theory_unc_hist(card_tool)

    np_map = {
            "Lambda2" : ["-0.25", "0.25",],
            "Lambda4" : [".01", ".16",],
            "Delta_Lambda2" : ["-0.02", "0.02",]
    }

    for nuisance,vals in np_map.items():
        card_tool.addSystematic(name=theory_hist,
            processes=samples,
            group="resumNonpert",
            systAxes=["vars",],
            passToFakes=to_fakes,
            action=lambda h,n=nuisance,vs=vals: h[{"vars" : [n+v for v in vs]}],
            outNames=[f"scetlibNP{nuisance}Down", f"scetlibNP{nuisance}Up"],
            rename=f"scetlib{nuisance}NP{name_append}",
        )

def theory_unc_hist(card_tool, unc=""):
    theory_unc = input_tools.args_from_metadata(card_tool, "theoryCorr")
    if not theory_unc:
        logger.error("Can not add resummation uncertainties. No theory correction was applied!")
    if theory_unc[0] != "scetlib_dyturbo":
        raise ValueError(f"The theory uncertainty hist {theory_unc} doesn't have the resummation uncertainty implemented")

    return theory_unc[0]+("Corr" if not unc else unc)

def add_common_np_uncertainties(card_tool, samples, to_fakes):

    # NOTE: The map needs to be keyed on the base procs not the group names, which is
    # admittedly a bit nasty
    expanded_samples = card_tool.datagroups.getProcNames(samples)
    logger.debug(f"expanded_samples: {expanded_samples}")
    obs = card_tool.project[:]
    if not obs:
        raise ValueError("Failed to find the observable names for the resummation uncertainties")

    theory_hist = theory_unc_hist(card_tool)
    # TODO: For now these are separate hists, to rerun everything together in the future
    if "scetlib_dyturbo" not in theory_hist:
        raise ValueError("Did not find new scetlib non-perturbative correction")

    # Since "c_nu = 0.1 is the central value, it doesn't show up in the name"
    ranges = {"gamma" : ["c_nu-0.1-omega_nu0.5", "omega_nu0.5", ]}
    for np_nuisance, npvars in ranges.items():
        nuisance_name = f"scetlibNP{np_nuisance}"
        card_tool.addSystematic(name=theory_hist,
            processes=samples,
            passToFakes=to_fakes,
            systAxes=["vars"],
            action=lambda h,np=npvars: h[{"vars" : npvars}],
            outNames=[f"{nuisance_name}Down", f"{nuisance_name}Up"],
def add_pdf_uncertainty(card_tool, samples, to_fakes, action=None, from_corr=False, scale=1):
    pdf = input_tools.args_from_metadata(card_tool, "pdfs")[0]
    pdfInfo = theory_tools.pdf_info_map("ZmumuPostVFP", pdf)
    pdfName = pdfInfo["name"]
    pdf_hist = pdfName

    if from_corr:
        theory_unc = input_tools.args_from_metadata(card_tool, "theoryCorr")
        if not theory_unc:
            logger.error("Can not add resummation uncertainties. No theory correction was applied!")
        pdf_hist = f"scetlib_dyturbo{pdf.upper()}Vars" 
        if pdf_hist not in theory_unc:
            logger.error(f"Did not find {pdf_hist} correction in file! Cannot use SCETlib+DYTurbo PDF uncertainties")
        pdf_hist += "Corr"

    logger.info(f"Using PDF hist {pdf_hist}")

    pdf_ax = "vars" if from_corr else "pdfVar"
    symHessian = pdfInfo["combine"] == "symHessian"
    pdf_args = dict(
        processes=samples,
        mirror=True if symHessian else False,
        group=pdfName,
        passToFakes=to_fakes,
        actionMap=action,
        scale=pdfInfo.get("scale", 1)*scale,
        systAxes=[pdf_ax],
    )
    if from_corr:
        card_tool.addSystematic(pdf_hist, 
            outNames=[""]+theory_tools.pdfNamesAsymHessian(pdfInfo['entries'], pdfset=pdf.upper())[1:],
            **pdf_args
        )
    else:
        card_tool.addSystematic(pdf_hist, 
            skipEntries=[{pdf_ax : "^pdf0[a-z]*"}],
            **pdf_args
        )

    # TODO: For now only MiNNLO alpha_s is supported
    asRange = pdfInfo['alphasRange']
    card_tool.addSystematic(f"{pdfName}alphaS{asRange}", 
        processes=samples,
        mirror=False,
        group=pdfName,
        systAxes=["alphasVar"],
        systNameReplace=[("as", "pdfAlphaS")]+[("0116", "Down"), ("0120", "Up")] if asRange == "002" else [("0117", "Down"), ("0119", "Up")],
        scale=0.75, # TODO: this depends on the set, should be provided in theory_tools.py
        passToFakes=to_fakes,
    )


def add_recoil_uncertainty(card_tool, samples, passSystToFakes=False, pu_type="highPU", flavor="", group_compact=True):
    met = input_tools.args_from_metadata(card_tool, "met")
    if flavor == "":
        flavor = input_tools.args_from_metadata(card_tool, "flavor")
    rtag = f"{pu_type}_{flavor}_{met}"
    if not rtag in recoil_tools.recoil_cfg:
        logger.warning(f"Recoil corrections for {pu_type}, {flavor}, {met} not available.")
        return
    recoil_cfg = recoil_tools.recoil_cfg[rtag]
    recoil_vars = list(recoil_cfg['corr_z'].keys()) + list(recoil_cfg['unc_z'].keys())
    recoil_grps = recoil_vars
    if group_compact:
        recoil_grps = ["CMS_recoil"]*len(recoil_cfg)
    for i, tag in enumerate(recoil_vars):
        card_tool.addSystematic("recoilUnc_%s" % tag,
            processes=samples,
            mirror = False,
            group = recoil_grps[i],
            systAxes = ["recoilVar"],
            passToFakes=passSystToFakes,
        )
