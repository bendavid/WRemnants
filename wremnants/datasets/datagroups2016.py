from utilities import boostHistHelpers as hh, logging
from wremnants import histselections as sel


logger = logging.child_logger(__name__)
    
def make_datagroups_2016(dg, combine=False, pseudodata_pdfset = None, applySelection=True, excludeGroups=None, filterGroups=None, simultaneousABCD=False):
    # reset datagroups
    dg.groups = {}

    if dg.mode == "wmass":
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
    dg.addGroup("Zmumu",
        members = dg.get_members_from_results(startswith=["Zmumu"]),
        selectOp = sigOp,
    ) 
    dg.addGroup("Ztautau",
        members = dg.get_members_from_results(startswith=["Ztautau"]),
        selectOp = sigOp,
    )
    dg.addGroup("PhotonInduced",
        members = dg.get_members_from_results(startswith=["GG", "QG"]),
        selectOp = sigOp,
    )

    if pseudodata_pdfset and dg.combine:
        dg.addGroup(f"pdf{pseudodata_pdfset.upper()}_sum",
            label = f"pdf{pseudodata_pdfset.upper()}",
            color = "dimgray"
        )
    if dg.mode in ["vgen", "wmass"]:
        dg.addGroup("Wmunu",
            members = dg.get_members_from_results(startswith=["Wplusmunu", "Wminusmunu"]),
            selectOp = sigOp,
        )
        dg.addGroup("Wtaunu",
            members = dg.get_members_from_results(startswith=["Wplustaunu", "Wminustaunu"]),
            selectOp = sigOp,
        )
        dg.addGroup("DYlowMass",
            members = dg.get_members_from_results(startswith=["DYlowMass", "DYJetsToMuMuMass10to50"]),
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
        dg.addGroup("QCD",
            members = dg.get_members_from_results(startswith=["QCD"]),
            selectOp = sigOp,
        )   
    else:
        dg.addGroup("Other",
            members = dg.get_members_from_results(not_startswith=["Zmumu", "Ztautau", "QCD", "GG", "QG"]),
        )

    dg.filterGroups(filterGroups)
    dg.excludeGroups(excludeGroups)

    if dg.mode == "wmass":
        # add all processes to the fake contributions after filtered and excluded groups
        dg.addGroup("Fake",
            members = [member for sublist in [v.members for k, v in dg.groups.items() if k != "QCD"] for member in sublist],
            scale = lambda x: 1. if x.is_data else -1,
            selectOp = fakeOp,
            selectOpArgs = fakeOpArgs
        )
        dg.filterGroups(filterGroups)
        dg.excludeGroups(excludeGroups)

    return dg
