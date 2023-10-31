from utilities import logging,common, boostHistHelpers as hh
from utilities.io_tools import input_tools
from wremnants import syst_tools,theory_tools
import numpy as np
import re

logger = logging.child_logger(__name__)

class TheoryHelper(object):
    valid_np_models = ["Lambda", "Omega", "Delta_Lambda", "Delta_Omega", "binned_Omega", "none"]
    def __init__(self, card_tool):
        for group in ['signal_samples', 'signal_samples_inctau', 'single_v_samples']:
            if group not in card_tool.procGroups:
                raise ValueError(f"Must define '{group}' procGroup in CardTool for theory uncertainties")
        
        self.card_tool = card_tool
        corr_hists = input_tools.args_from_metadata(self.card_tool, "theoryCorr")
        self.corr_hist_name = (corr_hists[0]+"Corr") if corr_hists else None
        self.syst_ax = "vars"
        self.corr_hist = None
        self.resumUnc = None
        self.np_model = "Delta_Lambda"
        self.pdf_from_corr = False
        self.scale_pdf_unc = 1.
        self.tnp_magnitude = 1.
        self.mirror_tnp = True
        self.minnlo_unc = 'byHelicityPt'

    def sample_label(self, sample_group):
        if sample_group not in self.card_tool.procGroups:
            raise ValueError(f"Failed to find sample group {sample_group} in predefined groups")
        if not self.card_tool.procGroups[sample_group]:
            logger.warning(f"Sample group {sample_group} is empty")

        return self.card_tool.procGroups[sample_group][0][0] 

    def configure(self, resumUnc, np_model,
            propagate_to_fakes=True, 
            tnp_magnitude=1,
            tnp_scale=1.,
            mirror_tnp=True,
            pdf_from_corr=False,
            pdf_action=None,
            scale_pdf_unc=1.,
            minnlo_unc='byHelicityPt'):

        self.set_resum_unc_type(resumUnc)
        self.set_np_model(np_model)
        self.set_propagate_to_fakes(propagate_to_fakes)
        self.set_minnlo_unc(minnlo_unc)

        self.tnp_magnitude = tnp_magnitude
        self.tnp_scale = tnp_scale
        self.mirror_tnp = mirror_tnp
        self.pdf_from_corr = pdf_from_corr
        self.pdf_action = pdf_action
        self.scale_pdf_unc = scale_pdf_unc
        self.minnlo_unc = minnlo_unc
        self.samples = []

    def add_all_theory_unc(self, nonsig=True):
        self.samples = ["signal_samples_inctau", "single_v_nonsig_samples"] if nonsig else ["signal_samples"]
        self.add_nonpert_unc(model=self.np_model)
        self.add_resum_unc(magnitude=self.tnp_magnitude, mirror=self.mirror_tnp, scale=self.tnp_scale)
        self.add_pdf_uncertainty(from_corr=self.pdf_from_corr, action=self.pdf_action, scale=self.scale_pdf_unc)

    def set_minnlo_unc(self, minnloUnc):
        self.minnlo_unc = minnloUnc

    def set_resum_unc_type(self, resumUnc):
        if not resumUnc or resumUnc == "none":
            self.resumUnc = None
            return

        if not self.corr_hist_name:
            raise ValueError("Cannot add resummation uncertainties. No theory correction was applied!")

        if input_tools.args_from_metadata(self.card_tool, "theoryCorrAltOnly"):
            raise ValueError("The theory correction was only applied as an alternate hist. Using it for systs isn't well defined!")

        signal_samples = self.card_tool.procGroups['signal_samples']
        self.corr_hist = self.card_tool.getHistsForProcAndSyst(signal_samples[0], self.corr_hist_name)

        if resumUnc == "tnp":
            self.tnp_nuisances = self.card_tool.match_str_axis_entries(self.corr_hist.axes[self.syst_ax], 
                                    ["^gamma_.*[+|-]\d+", "^b_.*[+|-]\d+", "^s[+|-]\d+", "^h_.*\d+"])
            if not self.tnp_nuisances:
                raise ValueError(f"Did not find TNP uncertainties in hist {self.corr_hist_name}")

            self.tnp_names = set([x.split("+")[0].split("-")[0] for x in self.tnp_nuisances])
            
        self.resumUnc = resumUnc
        
    def add_resum_unc(self, magnitude=1, mirror=False, scale=1):
        if not self.resumUnc:
            logger.warning("No resummation uncertainty will be applied!")

        if self.resumUnc == "tnp":
            self.add_resum_tnp_unc(magnitude, mirror, scale)

        if self.minnlo_unc and self.minnlo_unc not in ["none", None]:
            for sample_group in self.samples:
                if self.card_tool.procGroups.get(sample_group, None):
                    self.add_minnlo_scale_uncertainty(sample_group, rebin_pt=common.ptV_binning[::2])

    def add_minnlo_scale_uncertainty(self, sample_group, use_hel_hist=False, rebin_pt=None):
        if not sample_group or sample_group not in self.card_tool.procGroups:
            logger.warning(f"Skipping QCD scale syst '{self.minnlo_unc}' for group '{sample_group}.' No process to apply it to")
            return
            
        helicity = "Helicity" in self.minnlo_unc
        pt_binned = "Pt" in self.minnlo_unc
        scale_hist = "qcdScale" if not (helicity or use_hel_hist) else "qcdScaleByHelicity"
        if "helicity" in scale_hist.lower():
            use_hel_hist = True

        # All possible syst_axes
        # TODO: Move the axes to common and refer to axis_chargeVgen etc by their name attribute, not just
        # assuming the name is unchanged
        obs = self.card_tool.fit_axes
        pt_ax = "ptVgen" if "ptVgen" not in obs else "ptVgenAlt"
        charge_ax = "chargeVgen" if "chargeVgen" not in obs else "chargeVgenAlt"
        sum_axes = [pt_ax, charge_ax]

        syst_axes = sum_axes + ["muRfact", "muFfact"]
        syst_ax_labels = ["PtV", "genQ", "muR", "muF"]
        format_with_values = ["low", "center", "center", "center"]
        if use_hel_hist:
            syst_axes.insert(2, "helicity")
            syst_ax_labels.insert(2, "AngCoeff")
            format_with_values.insert(2, "low")

        group_name = f"QCDscale{self.sample_label(sample_group)}"
        # Exclude all combinations where muR = muF = 1 (nominal) or where
        # they are extreme values (ratio = 4 or 1/4)
        skip_entries = [{"muRfact" : 1.j, "muFfact" :  1.j}, {"muRfact" : 0.5j, "muFfact" : 2.j}, 
                        {"muRfact" : 2.j, "muFfact" : 0.5j}]
        # In order to make prettier names than the automated ones.
        # No harm in leaving extra replaces that won't be triggered
        name_replace = [("muR2muF2", "muRmuFUp"), ("muR0p5muF0p5", "muRmuFDown"), ("muR2muF1", "muRUp"), 
                            ("muR0p5muF1", "muRDown"), ("muR1muF0p5", "muFDown"), ("muR1muF2", "muFUp"),
        ]
        action_map = {}

        # NOTE: The map needs to be keyed on the base procs not the group names, which is
        # admittedly a bit nasty
        expanded_samples = self.card_tool.getProcNames([sample_group])
        logger.debug(f"using {scale_hist} histogram for QCD scale systematics")
        logger.debug(f"expanded_samples: {expanded_samples}")
        func = syst_tools.scale_helicity_hist_to_variations 
        action_args = {"sum_axes" : sum_axes}
        if self.card_tool.datagroups.mode == "vgen":
            func = syst_tools.gen_scale_helicity_hist_to_variations
            action_args["gen_obs"] = obs
            action_args["sum_axes"].extend(["y", "massVgen"])
            action_args["pt_ax"] = pt_ax
            action_args["gen_axes"] = syst_axes[:-2] #TODO: Might want to incorporate support for yVgen also

        action_map = {proc : func for proc in expanded_samples}
            
        # Determine if it should be summed over based on minnlo_unc. If not,
        # Remove it from the sum list and set names appropriately
        def set_sum_over_axis(identifier, ax_name):
            nonlocal sum_axes,group_name,syst_axes,syst_ax_labels
            if identifier in self.minnlo_unc:
                sum_axes.remove(ax_name)
                group_name += identifier
            elif ax_name in syst_axes:
                idx = syst_axes.index(ax_name)
                syst_axes.pop(idx)
                syst_ax_labels.pop(idx)
                format_with_values.pop(idx)

        for ax,name in zip(sum_axes[:], ["Pt", "Charge", "Helicity"]):
            set_sum_over_axis(name, ax)

        if pt_binned:
            signal_samples = self.card_tool.procGroups['signal_samples']
            binning = np.array(rebin_pt) if rebin_pt else None

            hscale = self.card_tool.getHistsForProcAndSyst(signal_samples[0], scale_hist)
            # A bit janky, but refer to the original ptVgen ax since the alt hasn't been added yet
            orig_binning = hscale.axes[pt_ax.replace("Alt", "")].edges
            if not hh.compatibleBins(orig_binning, binning):
                logger.warning(f"Requested binning {binning} is not compatible with hist binning {orig_binning}. Will not rebin!")
                binning = orig_binning

            if self.resumUnc:
                pt30_idx = np.argmax(binning > 30)

                if helicity:
                    # Drop the uncertainties for < 30 for sigma_-1
                    skip_entries.extend([{"helicity" : -1.j, pt_ax : complex(0, x)} for x in binning[:pt30_idx-1]])
                else:
                    # Drop the uncertainties for < 30
                    skip_entries.extend([{pt_ax : complex(0, x)} for x in binning[:pt30_idx-1]])

            action_args["rebinPtV"] = binning

        if helicity:
            # Drop the uncertainty of A5,A6,A7
            skip_entries.extend([{"helicity" : complex(0, i)} for i in (5,6,7)])

        # Skip MiNNLO unc. 
        if self.resumUnc and not (pt_binned or helicity):
            logger.warning("Without pT or helicity splitting, only the SCETlib uncertainty will be applied!")
        else:
            #FIXME Maybe put W and Z nuisances in the same group
            group_name += f"MiNNLO"
            self.card_tool.addSystematic(scale_hist,
                actionMap=action_map,
                actionArgs=action_args,
                processes=[sample_group],
                group=group_name,
                splitGroup={"QCDscale": ".*"},
                systAxes=syst_axes,
                labelsByAxis=syst_ax_labels,
                # Exclude all combinations where muR = muF = 1 (nominal) or where
                # they are extreme values (ratio = 4 or 1/4)
                skipEntries=skip_entries,
                systNameReplace=name_replace,
                baseName=group_name+"_",
                formatWithValue=format_with_values,
                passToFakes=self.propagate_to_fakes,
                rename=group_name, # Needed to allow it to be called multiple times
            )

    def set_propagate_to_fakes(self, to_fakes):
        self.propagate_to_fakes = to_fakes

    def add_resum_tnp_unc(self, magnitude, mirror, scale=1):
        syst_ax = self.corr_hist.axes[self.syst_ax]

        np_mag = lambda x: x.split("+")[-1].split("-")[-1]
        tnp_magnitudes = set([np_mag(x) for x in self.tnp_nuisances])
        tnp_has_mirror = ("+" in self.tnp_nuisances[0] and self.tnp_nuisances[0].replace("+", "-") in syst_ax) or \
                            ("-" in self.tnp_nuisances[0] and self.tnp_nuisances[0].replace("-", "+") in syst_ax)

        if not any(np.isclose(float(x), magnitude) for x in tnp_magnitudes):
            raise ValueError(f"TNP magnitude variation {magnitude} is not present in the histogram {self.corr_hist_name}. Options are {tnp_magnitudes}")

        qqV_mag = magnitude/2.
        if qqV_mag == 2.5:
            qqV_mag = 2.0

        selected_tnp_nuisances = [x for x in self.tnp_nuisances if np.isclose(magnitude, float(np_mag(x))) or ("h_qqV" in x and np.isclose(qqV_mag, float(np_mag(x))))]
        central_var = syst_ax[0]
        
        if mirror:
            name_replace = []
            if tnp_has_mirror:
                selected_tnp_nuisances = [x for x in selected_tnp_nuisances if "+" in x]
        else:
            if not tnp_has_mirror:
                raise ValueError(f"Uncertainty hist {self.corr_hist_name} does not have double sided TNP nuisances. " \
                                "mirror=False is therefore not defined")
            name_replace = [(f"+{m}", "Up") for m in tnp_magnitudes]+[(f"-{m}", "Up") for m in tnp_magnitudes]

        logger.debug(f"TNP nuisances in correction hist: {self.tnp_nuisances}")
        logger.debug(f"Selected TNP nuisances: {selected_tnp_nuisances}")

        self.card_tool.addSystematic(name=self.corr_hist_name,
            processes=['single_v_samples'],
            group="resumTNP",
            splitGroup={"resum": ".*"},
            systAxes=["vars"],
            passToFakes=self.propagate_to_fakes,
            systNameReplace=name_replace,
            action=lambda h: h[{self.syst_ax : [central_var, *selected_tnp_nuisances]}],
            doActionBeforeMirror=True,
            mirror=mirror,
            scale=scale,
            skipEntries=[{self.syst_ax : central_var},],
            rename=f"resumTNP",
            systNamePrepend=f"resumTNP_",
        )

    def add_nonpert_unc(self, model):
        if not self.resumUnc:
            return

        self.set_np_model(model)
        if self.np_model:
            self.add_gamma_np_uncertainties()
            self.add_uncorrelated_np_uncertainties()
        else:
            logger.warning("Will not add any nonperturbative uncertainty!")

    def set_np_model(self, model):
        if model in ["none", None]:
            self.np_model = None
            return

        if model not in TheoryHelper.valid_np_models:
            raise ValueError(f"Model choice {model} is not a supported model. Valid choices are {TheoryHelper.valid_np_models}")
    
        signal_samples = self.card_tool.procGroups['signal_samples']
        if model == "binned_Omega":
            self.np_hist_name = self.corr_hist_name.replace("Corr", "Omega")
            self.np_hist = self.card_tool.getHistsForProcAndSyst(signal_samples[0], self.np_hist_name)
        else:
            self.np_hist_name = self.corr_hist_name
            self.np_hist = self.corr_hist if self.corr_hist else self.card_tool.getHistsForProcAndSyst(signal_samples[0], self.corr_hist_name)

        var_name = model.replace("binned_", "")

        if not any(var_name in x for x in self.np_hist.axes[self.syst_ax]):
            raise ValueError(f"NP model choice was '{model}' but did not find corresponding variations in the histogram")

        self.np_model = model

    def add_gamma_np_uncertainties(self):
        # Since "c_nu = 0.1 is the central value, it doesn't show up in the name"
        gamma_vals = list(filter(lambda x: x in self.corr_hist.axes[self.syst_ax], 
            ["c_nu-0.1-omega_nu0.5", "omega_nu0.5", "c_nu-0.25", "c_nu-0.25"]))

        if len(gamma_vals) != 2:
            raise ValueError(f"Failed to find consistent variation for gamma NP in hist {self.corr_hist_name}")
        logger.debug(f"Adding gamma uncertainties from syst entries {gamma_vals}")

        nuisance_name = f"scetlibNPgamma"
        self.card_tool.addSystematic(name=self.corr_hist_name,
            processes=["single_v_samples"],
            passToFakes=self.propagate_to_fakes,
            systAxes=[self.syst_ax],
            action=lambda h: h[{self.syst_ax : gamma_vals}],
            outNames=[f"{nuisance_name}Down", f"{nuisance_name}Up"],
            group="resumNonpert",
            splitGroup={"resum": ".*"},
            rename=nuisance_name,
        )

    def add_resum_scale_uncertainty():
        obs = self.card_tool.fit_axes[:]
        if not obs:
            raise ValueError("Failed to find the observable names for the resummation uncertainties")
        
        theory_hist = self.card_tool.getHistsForProcAndSyst(samples[0], theory_hist_name)
        resumscale_nuisances = match_str_axis_entries(h.axes[syst_ax], ["^nuB.*", "nuS.*", "^muB.*", "^muS.*",])

        expanded_samples = card_tool.datagroups.getProcNames(samples)
        syst_ax = "vars"

        card_tool.addSystematic(name=theory_hist,
            processes=samples,
            group="resumScale",
            splitGroup={"resum": ".*"},
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
            splitGroup={"resum": ".*"},
            passToFakes=to_fakes,
            systAxes=["vars"],
            actionMap={s : lambda h: h[{"vars" : ["kappaFO0.5-kappaf2.", "kappaFO2.-kappaf0.5", "mufdown", "mufup",]}] for s in expanded_samples},
            outNames=[f"scetlib_kappa{name_append}Up", f"scetlib_kappa{name_append}Down", f"scetlib_muF{name_append}Up", f"scetlib_muF{name_append}Down"],
            rename=f"resumFOScale{name_append}",
            systNamePrepend=f"resumScale{name_append}_",
        )

    def add_uncorrelated_np_uncertainties(self):
        np_map = {
            "Lambda2" : ["-0.25", "0.25",],
            "Lambda4" : [".01", ".16",],
            "Delta_Lambda2" : ["-0.02", "0.02",]
        } if "Lambda" in self.np_model else {
            "Omega" : ["0.", "0.8"],
            "Delta_Omega" : ["-0.02", "0.02"],
        }

        if "Delta" not in self.resumUnc:
            to_remove = list(filter(lambda x: "Delta" in x, np_map.keys()))
            for k in to_remove:
                np_map.pop(k)

        central_var = self.np_hist.axes[self.syst_ax][0]

        for label,vals in np_map.items():
            if not all(label+v in self.np_hist.axes[self.syst_ax] for v in vals):
                tmpvals = [x.replace(label, "") for x in self.np_hist.axes[self.syst_ax] if re.match(f"^{label}\d+", x)]
                if tmpvals:
                    logger.warning(f"Using variations {tmpvals} rather than default values {vals}!")
                    np_map[label] = tmpvals
                else:
                    raise ValueError(f"Failed to find all vars {vals} for var {label} in hist {self.np_hist_name}")

        for sample_group in self.samples:
            if not self.card_tool.procGroups.get(sample_group, None):
                continue
            label = self.sample_label(sample_group)
            for nuisance,vals in np_map.items():
                entries = [nuisance+v for v in vals]
                binned = "binned" in self.np_model
                if binned:
                    action=lambda h,e=entries: syst_tools.hist_to_variations(h[{self.syst_ax : [central_var, *e]}])
                else:
                    action=lambda h,e=entries: h[{self.syst_ax : e}]
                rename = f"scetlibNP{label}{nuisance}"
                self.card_tool.addSystematic(name=self.np_hist_name,
                    processes=[sample_group],
                    group="resumNonpert",
                    splitGroup={"resum": ".*"},
                    systAxes=[self.syst_ax] if not binned else ["absYVgenNP", "chargeVgenNP", "vars"],
                    passToFakes=self.propagate_to_fakes,
                    action=action,
                    outNames=[f"{rename}Down", f"{rename}Up"] if not binned else None,
                    systNameReplace=[(entries[1], f"{rename}Up"), (entries[0], f"{rename}Down"), ],
                    skipEntries=[{self.syst_ax : central_var}],
                    rename=rename,
                )

    def add_pdf_uncertainty(self, action=None, from_corr=False, scale=1):
        pdf = input_tools.args_from_metadata(self.card_tool, "pdfs")[0]
        pdfInfo = theory_tools.pdf_info_map("ZmumuPostVFP", pdf)
        pdfName = pdfInfo["name"]
        pdf_hist = pdfName

        if from_corr:
            theory_unc = input_tools.args_from_metadata(self.card_tool, "theoryCorr")
            if not theory_unc:
                logger.error("Can not add resummation uncertainties. No theory correction was applied!")
            pdf_hist = f"scetlib_dyturbo{pdf.upper()}Vars" 
            if pdf_hist not in theory_unc:
                logger.error(f"Did not find {pdf_hist} correction in file! Cannot use SCETlib+DYTurbo PDF uncertainties")
            pdf_hist += "Corr"

        logger.info(f"Using PDF hist {pdf_hist}")

        pdf_ax = self.syst_ax if from_corr else "pdfVar"
        symHessian = pdfInfo["combine"] == "symHessian"
        pdf_args = dict(
            processes=["single_v_samples"],
            mirror=True if symHessian else False,
            group=pdfName,
            splitGroup={f"{pdfName}NoAlphaS": '.*'},
            passToFakes=self.propagate_to_fakes,
            actionMap=action,
            scale=pdfInfo.get("scale", 1)*scale,
            systAxes=[pdf_ax],
        )
        if from_corr:
            self.card_tool.addSystematic(pdf_hist, 
                outNames=[""]+theory_tools.pdfNamesAsymHessian(pdfInfo['entries'], pdfset=pdf.upper())[1:],
                **pdf_args
            )
        else:
            self.card_tool.addSystematic(pdf_hist, 
                skipEntries=[{pdf_ax : "^pdf0[a-z]*"}],
                **pdf_args
            )

        # TODO: For now only MiNNLO alpha_s is supported
        asRange = pdfInfo['alphasRange']
        self.card_tool.addSystematic(f"{pdfName}alphaS{asRange}", 
            processes=["single_v_samples"],
            mirror=False,
            group=pdfName,
            splitGroup={f"{pdfName}AlphaS": '.*'},
            systAxes=["alphasVar"],
            skipEntries=[{"alphasVar" : "as0118"}],
            systNameReplace=[("as", "pdfAlphaS")]+[("0116", "Down"), ("0120", "Up")] if asRange == "002" else [("0117", "Down"), ("0119", "Up")],
            scale=0.75 if asRange == "002" else 1.5,
            passToFakes=self.propagate_to_fakes,
        )

    def add_resum_transition_uncertainty(self):
        obs = self.card_tool.fit_axes[:]

        for sample_group in self.samples:
            if self.card_tool.procGroups.get(sample_group, None):
                continue
            expanded_samples = self.card_tool.getProcNames([sample_group])
            name_append = self.sample_label(sample_group)

            self.card_tool.addSystematic(name=self.corr_hist_name,
                processes=[sample_group],
                group="resumTransition",
                splitGroup={"resum": ".*"},
                systAxes=["downUpVar"],
                passToFakes=to_fakes,
                # NOTE: I don't actually remember why this used no_flow=ptVgen previously, I don't think there's any harm in not using it...
                actionMap={s : lambda h: hh.syst_min_and_max_env_hist(h, obs, self.syst_ax, 
                    [x for x in h.axes["vars"] if "transition_point" in x]) for s in expanded_samples},
                outNames=[f"resumTransition{name_append}Up", f"resumTransition{name_append}Down"],
                rename=f"scetlibResumTransition{name_append}",
            )
