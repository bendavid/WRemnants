from utilities import common
from wremnants import syst_tools
import numpy as np
import logging

def add_scale_uncertainty(card_tool, scale_type, samples, to_fakes, pdf, name_append="", scetlib=False, use_hel_hist=False, rebin_pt=None):
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
        
        # TODO: Implement pT splitting for SCETlib
        nscetlib_vars=45
        common_args = dict(name=f"scetlib{pdf if pdf != 'nnpdf31' else ''}Corr_unc",
            processes=samples,
            group=group_name,
            systAxes=["vars"],
            passToFakes=to_fakes,
        )
        card_tool.addSystematic(**common_args,
            outNames=["", "resumFOScaleUp", "resumFOScaleDown"]+[""]*(nscetlib_vars-3),
            rename="scetlibCorr_resumFO", 
        )
        card_tool.addSystematic(**common_args,
            outNames=[""]*3+["resumLambdaDown", "resumLambdaUp"]+[""]*(nscetlib_vars-5),
            rename="scetlibCorr_resumLambda", 
        )

        common_args["systAxes"] = ["downUpVar"]
        common_args["labelsByAxis"] = ["downUpVar"]
        axes = ["recoil_reco"] if card_tool.histName == "reco_mll" else ["pt", "eta"] 
        card_tool.addSystematic(**common_args,
            # TODO: Should support other variables in the fit
            actionMap={s : lambda h: syst_tools.uncertainty_hist_from_envelope(h, axes, range(5, 9)) for s in expanded_samples},
            baseName="resumTransition", 
            rename="scetlibCorr_resumTrans", 
        )
        card_tool.addSystematic(**common_args,
            # TODO: Should support other variables in the fit
            actionMap={s : lambda h: syst_tools.uncertainty_hist_from_envelope(h, axes, range(9, 45)) for s in expanded_samples},
            baseName="resumScale", 
            rename="scetlibCorr_resumScale", 
        )

    if helicity:
        # Drop the uncertainty of A5,A6,A7
        skip_entries.extend([(*[-1]*(nsyst_dims-3), complex(i), -1, -1) for i in (5,6,7)])

    # Skip MiNNLO unc. 
    if scetlib and not (pt_binned or helicity):
        logging.warning("Without pT or helicity splitting, only the SCETlib uncertainty will be applied!")
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
