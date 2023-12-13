from utilities import boostHistHelpers as hh, logging
from wremnants import histselections as sel
import hist

logger = logging.child_logger(__name__)

def make_datagroups_lowPU(dg, combine=False, excludeGroups=None, filterGroups=None, applySelection=True, simultaneousABCD=False):
    # reset datagroups
    dg.groups = {}

    if dg.mode == "lowpu_w":
        fakeOpArgs = {"fakerate_integration_axes":[]}
        if applySelection:
            sigOp = sel.signalHistWmass
            fakeOp = sel.fakeHistABCD
        else:
            sigOp = None
            fakeOp = sel.fakeHistSimultaneousABCD
    else:
        sigOp = None
        fakeOp = None
        fakeOpArgs = None

    dg.addGroup("Data",
        members = dg.get_members_from_results(is_data=True),
        selectOp = sigOp,
    )
    dg.addGroup("Ztautau",
        members = dg.get_members_from_results(startswith="Ztautau"),
        selectOp = sigOp,
    )

    if dg.flavor == "mu" or dg.flavor == "mumu":  
        dg.addGroup("Zmumu",
            members = dg.get_members_from_results(startswith="Zmumu"),
            selectOp = sigOp,
        ) 
        if dg.mode == "lowpu_w":
            dg.addGroup("Wmunu",
                members = dg.get_members_from_results(startswith=["Wplusmunu", "Wminusmunu"]),
                selectOp = sigOp,
            )

    if dg.flavor == "e" or dg.flavor == "ee":  
        dg.addGroup("Zee",
            members = dg.get_members_from_results(startswith="Zee"),
            selectOp = sigOp,
        ) 
        if dg.mode == "lowpu_w":
            dg.addGroup("Wenu",
                members = dg.get_members_from_results(startswith=["Wplusenu", "Wminusenu"]),
                selectOp = sigOp,
            )

    if dg.mode == "lowpu_w":
        dg.addGroup("Wtaunu",
            members = dg.get_members_from_results(startswith=["Wplustaunu", "Wminustaunu"]),
            selectOp = sigOp,
        )
        dg.addGroup("Top",
            members = dg.get_members_from_results(startswith=["Top", "SingleT", "TT"]),
            selectOp = sigOp,
        )
        dg.addGroup("Diboson",
            members = dg.get_members_from_results(startswith=["Diboson", "WW", "WZ", "ZZ"]),
            selectOp = sigOp,
        )
    else:
        dg.addGroup("Other",
            members = dg.get_members_from_results(not_startswith=["Zmumu", "Zee", "Ztautau", "QCD"]),
        )

    dg.filterGroups(filterGroups)
    dg.excludeGroups(excludeGroups)

    if dg.mode == "lowpu_w":
        # add all processes to the fake contributions after filtered and excluded groups
        dg.addGroup(dg.fakeName,
            members = [member for sublist in [v.members for k, v in dg.groups.items() if k != "QCD"] for member in sublist],
            scale = lambda x: 1. if x.is_data else -1,
            selectOp = fakeOp,
            selectOpArgs = fakeOpArgs
        )
        dg.filterGroups(filterGroups)
        dg.excludeGroups(excludeGroups)


    return dg
