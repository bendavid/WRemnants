from utilities import boostHistHelpers as hh, common, logging, input_tools
from wremnants import syst_tools,theory_tools,recoil_tools
import numpy as np
import re

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
