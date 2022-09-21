from utilities import common
from wremnants import syst_tools

def add_scale_uncertainty(card_tool, scale_type, samples, to_fakes, name_append="", scetlib=False, use_hel_hist=False, rebin_pt=None):
    inclusiveScale = scale_type == "integrated"
    helicity = "Helicity" in scale_type

    scale_hist = "qcdScale" if not helicity else "qcdScaleByHelicity"
    # All possible syst_axes
    syst_axes = ["ptVgen", "chargeVgen", "muRfact", "muFfact"]
    syst_ax_labels = ["PtVBin", "genQ", "muR", "muF"]
    if helicity:
        syst_axes.insert(2, "Coeff")
        syst_ax_labels.insert(2, "helicity")

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
    if use_hel_hist:
        sum_axes.append("helicity")

    action_map = {proc : syst_tools.scale_helicity_hist_to_variations for proc in samples}
        
    # Determine if it should be summed over based on scale_type passed in. If not,
    # Remove it from the sum list and set names appropriately
    def set_sum_over_axis(identifier, ax_name, group_name_app):
        nonlocal sum_axes,group_name,syst_axes,syst_ax_labels
        if identifier in scale_type:
            sum_axes.remove(ax_name)
            group_name += group_name_app
        elif ax_name in syst_axes:
            idx = syst_axes.index(ax_name)
            syst_axes.pop(idx)
            syst_ax_labels.pop(idx)

    # TODO: Move the axes to common and refer to axis_chargeVgen etc by their name attribute, not just
    # assuming the name is unchanged
    for ax,iden,name in zip(sum_axes[:], ["Pt", "Charge", "Helicity"], ["ByPtV", "ByCharge", "ByHelicity"]):
        set_sum_over_axis(iden, ax, name)

    action_args = {"sum_axes" : sum_axes}
    if "Pt" in scale_type:
        action_args["rebinPtV"] = rebin_pt

    # Skip extreme muR/muF values for all bin combos (-1 = any)
    nsyst_dims = len(syst_axes)-len(sum_axes)
    skip_entries = [(*[-1]*nsyst_dims-2,*x) for x in skip_entries]

    if scetlib:
        # At least I hope this is the binning... Ideally would get it from the hist, but it's not trivial here
        binning = np.array(common.ptV_10quantiles_binning
        pt30_idx = np.argmax(binning > 30)
        # Drop the uncertainties for < 30
        if not helicity:
            skip_entries.extend([(complex(0,x),*[-1]*(nsyst_dim-1)) for x in binning[:pt30_idx]])
        else:
            nhel = 9
            other_dims = [([-1,i] if "chargeVgen" in syst_axes else [i])+[-1]*2 for i in range(nhel)]
            skip_entries.extend([(complex(0, x), *other) for other in other_dims for x in binning[:pt30_idx]])
        card_tool.addSystematic("scetlibCorr_unc",
            processes=samples,
            group="QCDscale_resum",
            systAxes=["vars"],
            systNames=["", "resumFOScaleUp", "resumFOScaleDown"]+[""]*41,
            passToFakes=to_fakes,
            rename="scetlibCorr_resumFO", 
        )
        card_tool.addSystematic("scetlibCorr_unc",
            processes=samples,
            group="QCDscale_lambda",
            systAxes=["vars"],
            systNames=[""]*3+["resumScaleLambdaDown", "resumScaleLambdaUp"]+[""]*39,
            passToFakes=to_fakes,
            rename="scetlibCorr_resumLambda", 
        )
        print("Here we are!")
        #card_tool.addSystematic("scetlibCorr_unc",
        #    processes=samples,
        #    group="QCDscale_lambda",
        #    systAxes=["vars"],
        #    systNames=[""]*5+["resumTransitionUp", "resumTransitionDown"]+[""]*37,
        #    passToFakes=to_fakes,
        #    rename="scetlibCorr_resumTrans", 
        #)

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
