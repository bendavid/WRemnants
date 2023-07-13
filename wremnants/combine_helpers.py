from utilities import boostHistHelpers as hh, common, logging, input_tools
from wremnants import syst_tools,theory_tools,recoil_tools
import numpy as np
import re

logger = logging.child_logger(__name__)

def add_modeling_uncertainty(card_tool, minnlo_scale, signal_samples, background_samples, to_fakes, resumType, wmass, scaleTNP=1, rebin_pt=None):
    scale_name = "W" if wmass else "Z"
    resum_samples = signal_samples+background_samples if wmass and resumType == "tnp" else signal_samples
    scale_label = scale_name if resumType != "tnp" else ""

    do_resum = resumType and resumType != "none"
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

def add_minnlo_scale_uncertainty(card_tool, scale_type, samples, to_fakes, resum, name_append="", use_hel_hist=False, rebin_pt=None):
    if not len(samples):
        logger.warning(f"Skipping QCD scale syst '{scale_type}', no process to apply it to")
        return
        
    helicity = "Helicity" in scale_type
    pt_binned = "Pt" in scale_type

    scale_hist = "qcdScale" if not (helicity or use_hel_hist) else "qcdScaleByHelicity"

    # All possible syst_axes
    # TODO: Move the axes to common and refer to axis_chargeVgen etc by their name attribute, not just
    # assuming the name is unchanged
    syst_axes = ["ptVgen", "chargeVgen", "muRfact", "muFfact"]
    syst_ax_labels = ["PtVBin", "genQ", "muR", "muF"]
    if helicity:
        syst_axes.insert(2, "helicity")
        syst_ax_labels.insert(2, "AngCoeff")

    group_name = f"QCDscale{name_append}"
    # Exclude all combinations where muR = muF = 1 (nominal) or where
    # they are extreme values (ratio = 4 or 1/4)
    skip_entries = [{"muRfact" : 1.j, "muFfact" :  1.j}, {"muRfact" : 0.5j, "muFfact" : 2.j}, 
                    {"muRfact" : 2.j, "muFfact" : 0.5j}]
    # In order to make prettier names than the automated ones.
    # No harm in leaving extra replaces that won't be triggered
    name_replace = [("muR2muF2", "muRmuFUp"), ("muR0muF0", "muRmuFDown"), ("muR2muF1", "muRUp"), 
                        ("muR0muF1", "muRDown"), ("muR1muF0", "muFDown"), ("muR1muF2", "muFUp"),
    ]
    action_map = {}
    sum_axes = ["ptVgen", "chargeVgen",]
    if use_hel_hist or helicity:
        sum_axes.append("helicity")

    # NOTE: The map needs to be keyed on the base procs not the group names, which is
    # admittedly a bit nasty
    expanded_samples = card_tool.datagroups.getProcNames(samples)
    logger.debug(f"using {scale_hist} histogram for QCD scale systematics")
    logger.debug(f"expanded_samples: {expanded_samples}")
    action_map = {proc : syst_tools.scale_helicity_hist_to_variations for proc in expanded_samples}
        
    # Determine if it should be summed over based on scale_type passed in. If not,
    # Remove it from the sum list and set names appropriately
    def set_sum_over_axis(identifier, ax_name):
        nonlocal sum_axes,group_name,syst_axes,syst_ax_labels
        if identifier in scale_type:
            sum_axes.remove(ax_name)
            group_name += identifier
        elif ax_name in syst_axes:
            idx = syst_axes.index(ax_name)
            syst_axes.pop(idx)
            syst_ax_labels.pop(idx)

    for ax,name in zip(sum_axes[:], ["Pt", "Charge", "Helicity"]):
        set_sum_over_axis(name, ax)

    action_args = {"sum_axes" : sum_axes}
    if pt_binned:
        action_args["rebinPtV"] = rebin_pt

    if helicity and resum:
        # Drop the uncertainty of A5,A6,A7
        skip_entries.extend([{"helicity" : complex(0, i)} for i in (5,6,7)])

    binning = np.array(common.ptV_10quantiles_binning)
    pt30_idx = np.argmax(binning > 30)
    if helicity:
        # Drop the uncertainties for < 30 for sigma_-1
        skip_entries.extend([{"helicity" : -1.j, "ptVgen" : complex(0, x)} for x in binning[:pt30_idx-1]])
    elif pt_binned:
        # Drop the uncertainties for < 30
        skip_entries.extend([{"ptVgen" : complex(0, x)} for x in binning[:pt30_idx-1]])

    # Skip MiNNLO unc. 
    if resum and not (pt_binned or helicity):
        logger.warning("Without pT or helicity splitting, only the SCETlib uncertainty will be applied!")
    else:
        #FIXME put W and Z nuisances in the same group
        group_name += f"MiNNLO"
        card_tool.addSystematic(scale_hist,
            actionMap=action_map,
            actionArgs=action_args,
            processes=samples,
            group=group_name,
            systAxes=syst_axes,
            labelsByAxis=syst_ax_labels,
            # Exclude all combinations where muR = muF = 1 (nominal) or where
            # they are extreme values (ratio = 4 or 1/4)
            skipEntries=skip_entries,
            systNameReplace=name_replace,
            baseName=group_name+"_",
            passToFakes=to_fakes,
            rename=group_name, # Needed to allow it to be called multiple times
        )


def add_resum_uncertainty(card_tool, samples, to_fakes, uncType, name_append="", scale=1):
    obs = card_tool.project[:]
    if not obs:
        raise ValueError("Failed to find the observable names for the resummation uncertainties")
    
    theory_hist = theory_unc_hist(card_tool)

    if input_tools.args_from_metadata(card_tool, "theoryCorrAltOnly"):
        logger.error("The theory correction was only applied as an alternate hist. Using its syst isn't well defined!")

    expanded_samples = card_tool.datagroups.getProcNames(samples)

    syst_ax = "vars"
    np_nuisances = ["^c_nu-*\d+", "^omega_nu-*\d+", "^Omega-*\d+"]
    both_exclude = ['^kappaFO.*','^recoil_scheme.*',"^transition_points.*",]+np_nuisances
    tnp_nuisance_names = ['gamma_cusp', 'gamma_mu_q', 'gamma_nu', 'h_qqV', 's', 'b_qqV', 'b_qqbarV', 'b_qqS', 'b_qqDS', 'b_qg', 'b_qg', ]
    tnp_nuisances = [x+"+5" if x != 'h_qqV' else x+"-2.0" for x in tnp_nuisance_names]
    resumscale_nuisances = ["^nuB.*", "nuS.*", "^muB.*", "^muS.*",]
    scale_nuisances = ["^mu.*", "^mu.*", "^nu",]
    if uncType == "tnp":
        # Exclude the scale uncertainty nuisances
        card_tool.addSystematic(name=theory_hist,
            processes=samples,
            group="resumTNP",
            systAxes=["vars"],
            passToFakes=to_fakes,
            systNameReplace=[("central", ""), ("pdf0", ""), ("+5", ""), ("-2.0", ""), ("+1", ""), ("+0.5", ""),],
            action=lambda h,tnp=tnp_nuisances: h[{syst_ax : ["central" if "central" in h.axes[syst_ax] else "pdf0"]+tnp_nuisances}],
                #("+1", "Up"), ("-1", "Down"), ("-0.5", "Down"), ("+0.5", "Up"), ("up", "Up"), ("down", "Down"),
                #("+1", "Up"), ("-1", "Down"), ("-0.5", "Down"), ("+0.5", "Up"), ("up", "Up"), ("down", "Down")],
            doActionBeforeMirror=True,
            mirror=True,
            scale=scale,
            skipEntries=[{syst_ax : "pdf0"},  {syst_ax : "central"}],
            rename=f"resumTNP{name_append}",
            systNamePrepend=f"resumTNP_",
        )
    else:
        # Exclude the tnp uncertainty nuisances
        card_tool.addSystematic(name=theory_hist,
            processes=samples,
            group="resumScale",
            passToFakes=to_fakes,
            skipEntries=[{syst_ax : x} for x in both_exclude+tnp_nuisances],
            systAxes=["downUpVar"], # Is added by the actionMap
            actionMap={s : lambda h: hh.syst_min_and_max_env_hist(h, obs, "vars", 
                [x for x in h.axes["vars"] if any(re.match(y, x) for y in resumscale_nuisances)]) for s in expanded_samples},
            outNames=[f"scetlibResumScale{name_append}Up", f"scetlibResumScale{name_append}Down"],
            rename=f"resumScale{name_append}",
            systNamePrepend=f"resumScale{name_append}_",
        )
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

<<<<<<< HEAD
def add_decorrelated_np_uncertainties(card_tool, samples, to_fakes, name_append, nuisances=["Omega"]):
    obs = card_tool.project[:]
    if not obs:
        raise ValueError("Failed to find the observable names for the resummation uncertainties")
=======
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
>>>>>>> b144363 (Fixes for new NP model)

    for nuisance,vals in np_map.items():
        card_tool.addSystematic(name=theory_hist,
            processes=samples,
            group="resumNonpert",
<<<<<<< HEAD
            systAxes=["absYVgenNP", "chargeVgenNP", "downUpVar"],
            passToFakes=to_fakes,
            actionMap={s : lambda h,np=np_nuisance: hh.syst_min_and_max_env_hist(syst_tools.hist_to_variations(h), obs, "vars",
                [x for x in h.axes["vars"] if re.match(f"^{np}-*\d+", x)]) for s in expanded_samples},
            baseName=f"scetlibNP{nuisance_name}_",
            rename=nuisance_name,
=======
            systAxes=["vars",],
            passToFakes=to_fakes,
            action=lambda h,n=nuisance,vs=vals: h[{"vars" : [n+v for v in vs]}],
            outNames=[f"scetlibNP{nuisance}Down", f"scetlibNP{nuisance}Up"],
            rename=f"scetlib{nuisance}NP{name_append}",
>>>>>>> b144363 (Fixes for new NP model)
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
<<<<<<< HEAD

    for np_nuisance in ["c_nu", "omega_nu"]:
        nuisance_name = f"scetlibNP{np_nuisance}"
        card_tool.addSystematic(name=theory_hist,
            processes=samples,
=======
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
>>>>>>> b144363 (Fixes for new NP model)
            group="resumNonpert",
            systAxes=["downUpVar"],
            passToFakes=to_fakes,
            actionMap={s : lambda h,np=np_nuisance: hh.syst_min_and_max_env_hist(h, obs, "vars",
                [x for x in h.axes["vars"] if re.match(f"^{np}-*\d+", x)]) for s in expanded_samples},
            outNames=[f"{nuisance_name}Up", f"{nuisance_name}Down"],
            rename=nuisance_name,
        )

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
