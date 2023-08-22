from utilities import boostHistHelpers as hh, common, logging, input_tools
from wremnants import syst_tools,theory_tools,recoil_tools, unfolding_tools

from wremnants import histselections as sel
import hist
import numpy as np
import uproot
import h5py

logger = logging.child_logger(__name__)

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

def setSimultaneousABCD(cardTool, variation_fakerate=0.5, variation_normalization_fake=0.1):
    # Having 1 process for fakes, for each bin 3 free floating parameters, 2 normalization for lowMT and highMT and one fakerate between iso and anti iso
    logger.info(f"Set processes for simultaneous ABCD fit")

    # expected fake contribution
    hist_fake = sum([group.hists[cardTool.nominalName] if name == cardTool.dataName else -1*group.hists[cardTool.nominalName] for name, group in cardTool.datagroups.groups.items()])

    # setting errors to 0
    hist_fake.view(flow=True)[...] = np.stack((hist_fake.values(flow=True), np.zeros_like(hist_fake.values(flow=True))), axis=-1)

    if common.passIsoName not in hist_fake.axes.name or common.passMTName not in hist_fake.axes.name:
        raise RuntimeError(f'{common.passIsoName} and {common.passMTName} expected to be found in histogram, but only have axes {hist_fake.axes.name}')

    # axes in the correct ordering
    axes = cardTool.project[:]
    axes += [ax for ax in [common.passIsoName, common.passMTName] if ax not in axes]

    if set(hist_fake.axes.name) != set(axes) or hist_fake.axes.name[-2] != common.passIsoName or hist_fake.axes.name[-1] != common.passMTName:
        logger.info(f"Axes in histogram '{hist_fake.axes.name}' are not the same as required '{axes}' or in a different order than expected, try to project")
        hist_fake = hist_fake.project(*axes)

    # set the expected values in the signal region
    hist_fake.values(flow=True)[...,1,1] = sel.fakeHistABCD(hist_fake).values(flow=True)
    
    fakename = cardTool.getFakeName()

    cardTool.datagroups.addGroup(fakename, label = "Nonprompt", color = "grey", members=[],) #TODO check if existing group can be used
    cardTool.datagroups.groups[fakename].hists[f"{cardTool.nominalName}"] = hist_fake

    bin_sizes = [ax.size for ax in hist_fake.axes if ax.name in axes and ax.name not in [common.passIsoName, common.passMTName]]
    axes = [ax.name for ax in hist_fake.axes if ax.name in axes and ax.name not in [common.passIsoName, common.passMTName]]

    for i in range(np.product(bin_sizes)):
        current_i = i
        ax_idx = []
        for num in reversed(bin_sizes):
            ax_idx.insert(0, current_i % num)
            current_i //= num

        bin_name = "_".join([f"{ax}{ax_idx[j]}" for j, ax in enumerate(axes)])

        logger.debug(f"Now at {i}/{np.product(bin_sizes)}: {bin_name}")

        other_indices = {ax: ax_idx[j] for j, ax in enumerate(axes)}

        n_failMT_failIso = hist_fake[{**common.failIso, **common.failMT, **other_indices}].value
        n_failMT_passIso = hist_fake[{**common.passIso, **common.failMT, **other_indices}].value
        n_passMT_failIso = hist_fake[{**common.failIso, **common.passMT, **other_indices}].value

        n_failMT = n_failMT_failIso + n_failMT_passIso
        fr = n_failMT_failIso / n_failMT
        n_passMT = n_passMT_failIso / fr

        # systematic variation for fake normalization
        for nameMT, mtCut, n in (("LowMT", common.failMT, n_failMT), ("HighMT", common.passMT, n_passMT)):

            hist_var = hist.Hist(*[*hist_fake.axes, common.down_up_axis], storage=hist.storage.Double())
            hist_var.view(flow=True)[...] = np.stack((hist_fake.values(flow=True), hist_fake.values(flow=True)), axis=-1)

            hist_var[{**common.failIso, **mtCut, **{"downUpVar":0}, **other_indices}] = (1-variation_normalization_fake) * n * fr
            hist_var[{**common.passIso, **mtCut, **{"downUpVar":0}, **other_indices}] = (1-variation_normalization_fake) * n * (1-fr)

            hist_var[{**common.failIso, **mtCut, **{"downUpVar":1}, **other_indices}] = (1+variation_normalization_fake) * n * fr
            hist_var[{**common.passIso, **mtCut, **{"downUpVar":1}, **other_indices}] = (1+variation_normalization_fake) * n * (1-fr)

            cardTool.datagroups.groups[fakename].hists[f"{cardTool.nominalName}_N{fakename}{nameMT}_{bin_name}"] = hist_var

            cardTool.addSystematic(f"N{fakename}{nameMT}_{bin_name}",
                processes=[fakename],
                group=f"{fakename}{nameMT}",
                noConstraint=True,
                outNames=[f"N{fakename}{nameMT}_{bin_name}Down", f"N{fakename}{nameMT}_{bin_name}Up"],
                systAxes=["downUpVar"],
                labelsByAxis=["downUpVar"],
            )

        # systematic variation for fakerate, should be smaller 1 and bigger 0
        diff = min(fr, 1-fr)
        frUp = fr + variation_fakerate * diff
        frDn = fr - variation_fakerate * diff

        hist_var = hist.Hist(*[*hist_fake.axes, common.down_up_axis], storage=hist.storage.Double())
        hist_var.view(flow=True)[...] = np.stack((hist_fake.values(flow=True), hist_fake.values(flow=True)), axis=-1)

        hist_var[{**common.failIso, **common.failMT, **{"downUpVar":0}, **other_indices}] = n_failMT * frDn
        hist_var[{**common.passIso, **common.failMT, **{"downUpVar":0}, **other_indices}] = n_failMT * (1-frDn)
        hist_var[{**common.failIso, **common.passMT, **{"downUpVar":0}, **other_indices}] = n_passMT * frDn
        hist_var[{**common.passIso, **common.passMT, **{"downUpVar":0}, **other_indices}] = n_passMT * (1-frDn)

        hist_var[{**common.failIso, **common.failMT, **{"downUpVar":1}, **other_indices}] = n_failMT * frUp
        hist_var[{**common.passIso, **common.failMT, **{"downUpVar":1}, **other_indices}] = n_failMT * (1-frUp)
        hist_var[{**common.failIso, **common.passMT, **{"downUpVar":1}, **other_indices}] = n_passMT * frUp
        hist_var[{**common.passIso, **common.passMT, **{"downUpVar":1}, **other_indices}] = n_passMT * (1-frUp)

        cardTool.datagroups.groups[fakename].hists[f"{cardTool.nominalName}_r{fakename}_{bin_name}"] = hist_var

        cardTool.addSystematic(f"r{fakename}_{bin_name}",
            processes=[fakename],
            group=f"r{fakename}",
            noConstraint=True,
            outNames=[f"r{fakename}_{bin_name}Down", f"r{fakename}_{bin_name}Up"],
            systAxes=["downUpVar"],
            labelsByAxis=["downUpVar"],
        )

def setTheoryFitData(cardTool, fitresult):
    logger.info(f"Prepare theory fit: load measured differential cross secction distribution and covariance matrix")

    poi_type="pmaskedexp"
    base_process = "W" if cardTool.datagroups.wmass else "Z"

    if fitresult.endswith(".root"):
        
        rfile = uproot.open(fitresult)
        df = unfolding_tools.get_results(rfile, poi_type)

        # write out unfolded data as 1D hist
        hist_xsec = hist.Hist(
            hist.axis.Regular(bins=len(df), start=0.5, stop=len(df)+0.5, underflow=False, overflow=False), storage=hist.storage.Weight())
        hist_xsec.view(flow=False)[...] = np.stack([df["value"].values, (df["err_total"].values)**2], axis=-1)

        cardTool.theoryFitData = hist_xsec.values(flow=False).flatten()

        # write out covariance as 2D hist
        cardTool.theoryFitDataCov = unfolding_tools.matrix_poi(rfile, poi_type, base_process=base_process, axes=cardTool.project).values(flow=False)
    elif fitresult.endswith(".hdf5"):
        hfile = h5py.File(fitresult, mode='r')

        cardTool.theoryFitData = hfile["pmaskedexp_outvals"][...]

        npoi = len(cardTool.theoryFitData)
        # make matrix between POIs only; assume POIs come first
        cardTool.theoryFitDataCov = hfile["pmaskedexp_outcov"][:npoi,:npoi]
    else:
        raise NotImplementedError(f"Unkown data type for fitresult {fitresult}")
    
    cardTool.theoryFit = True
