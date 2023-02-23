from utilities import common,boostHistHelpers as hh,common
from wremnants import syst_tools,theory_tools
import numpy as np
import logging
import re

logger = common.child_logger(__name__)

def add_scale_uncertainty(card_tool, scale_type, samples, to_fakes, name_append="", scetlib=None, use_hel_hist=False, rebin_pt=None):
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
    skip_entries = [(1.j, 1.j), (0.5j, 2.j), (2.j, 0.5j)]
    # In order to make prettier names than the automated ones.
    # No harm in leaving extra replaces that won't be triggered
    name_replace = [("muR2muF2", "muRmuFUp"), ("muR0muF0", "muRmuFDown"), ("muR2muF1", "muRUp"), 
                        ("muR0muF1", "muRDown"), ("muR1muF0", "muFDown"), ("muR1muF2", "muFUp"),
    ]
    scaleActionArgs = {}
    action_map = {}
    sum_axes = ["ptVgen", "chargeVgen",]
    if use_hel_hist or helicity:
        sum_axes.append("helicity")

    # NOTE: The map needs to be keyed on the base procs not the group names, which is
    # admittedly a bit nasty
    expanded_samples = card_tool.datagroups.getProcNames(samples)
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
    skip_entries = [(*[-1]*(nsyst_dims-2),*x) for x in skip_entries]

    if scetlib:
        group_name += "Resum"
        # At least I hope this is the binning... Ideally would get it from the hist, but it's not trivial here
        other_dims = []

        if pt_binned:
            binning = np.array(common.ptV_10quantiles_binning)
            pt30_idx = np.argmax(binning > 30)
            # Drop the uncertainties for < 30
            other_dims = [[-1]*(nsyst_dims-1)]
            if helicity:
                # Skip only the 1+cos^2 term
                other_dims = [(*[-1]*(nsyst_dims-4), -1.j, -1, -1)]

            skip_entries.extend([(complex(0, x), *other) for other in other_dims for x in binning[:pt30_idx-1]])
        elif helicity:
            skip_entries.append((*[-1]*(nsyst_dims-2), -1.j, -1, -1))
        
        obs = ("eta", "pt") if not card_tool.project else card_tool.project
        # TODO: Implement pT splitting for SCETlib
        # TODO: The hist used needs to be configurable

        theory_hist = args_from_metadata(card_tool.datagroups.getMetaInfo(), "theory_corr", "scetlib_dyturboCorr_unc")
        if "--theory_corr_alt_only" in card_tool.datagroups.getMetaInfo()["command"]:
            logger.warning("The theory correction was only applied as an alternate hist. Using its syst isn't well defined!")

        card_tool.addSystematic(name=f"scetlib_dyturboN4LLCorr_unc",
            processes=samples,
            group="scetlibTNP",
            systAxes=["vars"],
            passToFakes=to_fakes,
            systNameReplace=[("central", ""), ("+1", "Up"), ("-1", "Down"), ("-0.5", "Down"), ("+0.5", "Up"), ("up", "Up"), ("down", "Down")],
            skipEntries=[('^kappaFO.*',), ('^recoil_scheme.*',), ('^c_nu.*',), ('Omega\d.\d*',), ("^transition_points.*",)],
            rename=f"scetlibTNP", 
        )
        card_tool.addSystematic(name=f"scetlib_dyturboN4LLCorr_unc",
            processes=samples,
            group="scetlibNonpert",
            systAxes=["downUpVar"],
            passToFakes=to_fakes,
            actionMap={s : lambda h: hh.syst_min_and_max_env_hist(h, obs, "vars", 
                ['c_nu-0.15-omega_nu0.43', 'c_nu0.05', 'c_nu0.5-omega_nu0.15', 'c_nu-0.5-omega_nu0.37', 'Omega0.', 'Omega0.71',],
            no_flow=["ptVgen"]) for s in expanded_samples},
            outNames=["scetlibNonpertUp", "scetlibNonpertDown"],
            rename=f"scetlibNonpert", 
        )
        card_tool.addSystematic(name=f"scetlib_dyturboN4LLCorr_unc",
            processes=samples,
            group="resumTransition",
            systAxes=["downUpVar"],
            passToFakes=to_fakes,
            actionMap={s : lambda h: hh.syst_min_and_max_env_hist(h, obs, "vars", 
                ['transition_points0.4_0.75_1.1', 'transition_points0.2_0.45_0.7', 'transition_points0.4_0.55_0.7', 'transition_points0.2_0.65_1.1'],
            no_flow=["ptVgen"]) for s in expanded_samples},
            outNames=["resumTransitionUp", "resumTransitionDown"],
            rename=f"scetlibResumTransition", 
        )

    if helicity:
        # Drop the uncertainty of A5,A6,A7
        skip_entries.extend([(*[-1]*(nsyst_dims-3), complex(0, i), -1, -1) for i in (5,6,7)])

    # Skip MiNNLO unc. 
    if scetlib and not (pt_binned or helicity):
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

# TODO: It's a bit dangerous to not double check what the default really is
def args_from_metadata(meta_data, arg, default):
    if "command" not in meta_data.keys():
        raise ValueError(f"Failed to find command in meta_data (meta_data keys are {meta_data.keys()}")

    command = meta_data["command"]
    sys_args = np.array([x.strip("'") for x in command.split()])
    matching = (sys_args == f"--{arg}") | (sys_args == f"-{arg}")
    if not np.count_nonzero(matching):
        logger.warning(f"Did not find argument {arg}. Assuming the default value {default}")
        return [default]

    idx = np.argmax(matching)
    isflag = np.vectorize(lambda x: bool(re.match("^-+[a-z]", x)))(sys_args)
    isflag[:idx+1] = False
    # Select args until the next flag (this will break if you have a positional arg afterwards...)
    last_idx = np.argmax(isflag)
    if not last_idx:
        last_idx = len(sys_args)
    return sys_args[idx+1:last_idx]

def add_pdf_uncertainty(card_tool, samples, to_fakes):
    pdf = args_from_metadata(card_tool.datagroups.getMetaInfo(), "pdfs", "msht20")[0]
    logger.info(f"Using PDF {pdf}")
    pdfInfo = theory_tools.pdf_info_map("ZmumuPostVFP", pdf)
    pdfName = pdfInfo["name"]

    if pdfInfo["combine"] == "symHessian":
        card_tool.addSystematic(pdfName, 
            processes=samples,
            mirror=True,
            group=pdfName,
            systAxes=["pdfVar"],
            # Needs to be a tuple, since for multiple axis it would be (ax1, ax2, ax3)...
            # -1 means all possible values of the mirror axis
            skipEntries=[("^pdf0[a-z]*", -1)],
            passToFakes=to_fakes,
        )
    else:
        card_tool.addSystematic(pdfName, 
            processes=samples,
            mirror=False,
            group=pdfName,
            systAxes=["pdfVar"],
            skipEntries=[("^pdf0[a-z]*",)],
            passToFakes=to_fakes,
            scale=pdfInfo["scale"] if "scale" in pdfInfo else 1,
        )

    asRange = pdfInfo['alphasRange']
    card_tool.addSystematic(f"{pdfName}alphaS{asRange}", 
        processes=samples,
        mirror=False,
        group=pdfName,
        systAxes=["alphasVar"],
        systNameReplace=[("0116", "Down"), ("0120", "Up")] if asRange == "002" else [("0117", "Down"), ("0119", "Up")],
        scale=0.75, # TODO: this depends on the set, should be provided in theory_tools.py
        passToFakes=to_fakes,
    )
