from utilities import boostHistHelpers as hh, common, logging,input_tools
from wremnants import syst_tools,theory_tools
import numpy as np
import re

logger = logging.child_logger(__name__)

def add_scale_uncertainty(card_tool, scale_type, samples, to_fakes, name_append="", resum=None, use_hel_hist=False, rebin_pt=None):
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
    expanded_samples = card_tool.datagroups.getProcNames(samples, afterFilter=True)
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

    # Skip extreme muR/muF values for all bin combos (-1 = any)
    nsyst_dims = len(syst_axes)
    #skip_entries = [(*[-1]*(nsyst_dims-2),*x) for x in skip_entries]

    if resum:
        group_name += f"Resum{resum.capitalize()}"
        if pt_binned:
            binning = np.array(common.ptV_10quantiles_binning)
            pt30_idx = np.argmax(binning > 30)
            # Drop the uncertainties for < 30
            skip_entries.extend([{"ptVgen" : complex(0, x)} for x in binning[:pt30_idx-1]])
        if helicity:
            skip_entries.append({"helicity" : -1.j})
        
        obs = card_tool.project
        if not obs:
            raise ValueError("Failed to find the observable names for the resummation uncertainties")
        # TODO: Implement pT splitting for SCETlib
        
        theory_unc = input_tools.args_from_metadata(card_tool, "theoryCorr")
        if not theory_unc:
            logger.error("Can not add resummation uncertainties. No theory correction was applied!")
        theory_unc = theory_unc[0]+"Corr"
        if theory_unc != "scetlib_dyturboCorr":
            raise ValueError(f"The theory uncertainty hist {theory_unc} doesn't have the resummation {resum} uncertainty implemented")

        if input_tools.args_from_metadata(card_tool, "theoryCorrAltOnly"):
            logger.error("The theory correction was only applied as an alternate hist. Using its syst isn't well defined!")

        syst_ax = "vars"
        np_nuisances = ["^c_nu-*\d+", "^omega_nu-*\d+", "^Omega-*\d+"]
        both_exclude = ['^kappaFO.*','^recoil_scheme.*',"^transition_points.*",]+np_nuisances
        tnp_nuisances = ["^gamma_.*", "b_.*", "s+*", "s-*",]
        resumscale_nuisances = ["^nuB.*", "nuS.*", "^muB.*", "^muS.*",]
        scale_nuisances = ["^mu.*", "^mu.*", "^nu",]
        if resum == "tnp":
            # Exclude the scale uncertainty nuisances
            card_tool.addSystematic(name=theory_unc,
                processes=samples,
                group="resumTNP",
                systAxes=["vars"],
                passToFakes=to_fakes,
                systNameReplace=[("central", ""), ("pdf0", ""), ("+1", "Up"), ("-1", "Down"), ("-0.5", "Down"), ("+0.5", "Up"), ("up", "Up"), ("down", "Down")],
                skipEntries=[{syst_ax : x} for x in both_exclude+scale_nuisances],
                rename=f"resumTNP", 
            )
        else:
            # Exclude the tnp uncertainty nuisances
            card_tool.addSystematic(name=theory_unc,
                processes=samples,
                group="resumTNP",
                passToFakes=to_fakes,
                skipEntries=[{syst_ax : x} for x in both_exclude+tnp_nuisances],
                systAxes=["downUpVar"], # Is added by the actionMap
                actionMap={s : lambda h: hh.syst_min_and_max_env_hist(h, obs, "vars", 
                    [x for x in h.axes["vars"] if any(re.match(y, x) for y in resumscale_nuisances)]) for s in expanded_samples},
                outNames=["scetlibResumScaleUp", "scetlibResumScaleDown"],
                rename=f"resumScale", 
            )
            card_tool.addSystematic(name=theory_unc,
                processes=samples,
                group="resumScale",
                passToFakes=to_fakes,
                systAxes=["vars"],
                actionMap={s : lambda h: h[{"vars" : ["kappaFO0.5-kappaf2.", "kappaFO2.-kappaf0.5", "mufdown", "mufup",]}] for s in expanded_samples},
                outNames=["scetlib_kappaUp", "scetlib_kappaDown", "scetlib_muFUp", "scetlib_muFDown"],
                rename=f"resumFOScale", 
            )

        for np_nuisance in ["c_nu", "omega_nu", "Omega"]:
            nuisance_name = f"scetlibNP{np_nuisance}"
            card_tool.addSystematic(name=theory_unc,
                processes=samples,
                group="resumNonpert",
                systAxes=["downUpVar"],
                passToFakes=to_fakes,
                actionMap={s : lambda h: hh.syst_min_and_max_env_hist(h, obs, "vars", 
                    [x for x in h.axes["vars"] if re.match(f"^{np_nuisance}-*\d+", x)]) for s in expanded_samples},
                outNames=[f"{nuisance_name}Up", f"{nuisance_name}Down"],
                rename=nuisance_name, 
            )
        card_tool.addSystematic(name=theory_unc,
            processes=samples,
            group="resumTransition",
            systAxes=["downUpVar"],
            passToFakes=to_fakes,
            # NOTE: I don't actually remember why this used no_flow=ptVgen previously, I don't think there's any harm in not using it...
            actionMap={s : lambda h: hh.syst_min_and_max_env_hist(h, obs, "vars", 
                [x for x in h.axes["vars"] if "transition_point" in x]) for s in expanded_samples},
            outNames=["resumTransitionUp", "resumTransitionDown"],
            rename=f"scetlibResumTransition", 
        )

    if helicity:
        # Drop the uncertainty of A5,A6,A7
        skip_entries.extend([{"helicity" : complex(0, i)} for i in (5,6,7)])

    # Skip MiNNLO unc. 
    if resum and not (pt_binned or helicity):
        logger.warning("Without pT or helicity splitting, only the SCETlib uncertainty will be applied!")
    else:
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

def add_pdf_uncertainty(card_tool, samples, to_fakes, action=None):
    pdf = input_tools.args_from_metadata(card_tool, "pdfs")[0]
    logger.info(f"Using PDF {pdf}")
    pdfInfo = theory_tools.pdf_info_map("ZmumuPostVFP", pdf)
    pdfName = pdfInfo["name"]

    pdf_ax = "pdfVar"
    symHessian = pdfInfo["combine"] == "symHessian"
    card_tool.addSystematic(pdfName, 
        processes=samples,
        mirror=True if symHessian else False,
        group=pdfName,
        systAxes=[pdf_ax],
        # Needs to be a tuple, since for multiple axis it would be (ax1, ax2, ax3)...
        # -1 means all possible values of the mirror axis
        skipEntries=[{pdf_ax : "^pdf0[a-z]*"}],
        passToFakes=to_fakes,
        actionMap=action,
        scale=pdfInfo["scale"] if "scale" in pdfInfo else 1,
    )

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
