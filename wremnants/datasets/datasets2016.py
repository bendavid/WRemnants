import narf
import logging
import subprocess
import glob
import pathlib
import socket
import logging
#set the debug level for logging incase of full printout 
from wremnants.datasets.datasetDict_v9 import dataDictV9
from wremnants.datasets.datasetDict_v8 import *
from wremnants.datasets.datasetDict_gen import genDataDict

logger = logging.getLogger("wremnants").getChild(__name__.split(".")[-1])

lumicsv = f"{pathlib.Path(__file__).parent.parent}/data/bylsoutput.csv"
lumijson = f"{pathlib.Path(__file__).parent.parent}/data/Cert_271036-284044_13TeV_Legacy2016_Collisions16_JSON.txt"

#TODO add the rest of the samples!
def makeFilelist(paths, maxFiles=-1, format_args={}):
    filelist = []
    for path in paths:
        if format_args:
            path = path.format(**format_args)
            logger.debug(f"Reading files from path {path}")
        filelist.extend(glob.glob(path) if path[:4] != "/eos" else buildXrdFileList(path, "eoscms.cern.ch"))
    return filelist if maxFiles < 0 else filelist[:maxFiles]

## TODO: use consistent names for filt and excludeGroup, filt was implemented before, should call it filterGroup
def getDatasets(maxFiles=-1, filt=None, excludeGroup=None, mode=None, base_path=None, nanoVersion="v9", 
        data_tag="TrackFitV718_NanoProdv1", mc_tag="TrackFitV718_NanoProdv1"):
    if not base_path:
        hostname = socket.gethostname()
        if hostname == "lxplus8s10.cern.ch":
            base_path = "/scratch/shared/NanoAOD"
        elif "mit.edu" in hostname:
            base_path = "/scratch/submit/cms/wmass/NanoAOD"
        elif hostname == "cmsanalysis.pi.infn.it":
            base_path = "/scratchnvme/wmass/NANOV9/newNTuples" #temporary

    logger.info(f"Loading samples from {base_path}.")

    if nanoVersion == "v8":
        dataDict = dataDictV8
        logger.info('Using NanoAOD V8')
    elif nanoVersion == "v9":
        dataDict = dataDictV9
    else:
        raise ValueError("Only NanoAODv8 and NanoAODv9 are supported")

    if mode == "gen":
        dataDict.update(genDataDict)

    narf_datasets = []
    for sample,info in dataDict.items():
        if sample in genDataDict:
            base_path = base_path.replace("NanoAOD", "NanoGen")

        is_data = "data" in sample[:4]

        prod_tag = data_tag if is_data else mc_tag 
        paths = makeFilelist(info["filepaths"], maxFiles, format_args=dict(BASE_PATH=base_path, NANO_PROD_TAG=prod_tag))

        if not paths:
            logger.warning(f"Failed to find any files for dataset {sample}. Looking at {info['filepaths']}. Skipping!")
            continue

        narf_info = dict(
            name=sample,
            filepaths=paths,
        )

        if is_data:
            if mode == "gen":
                continue
            narf_info.update(dict(
                is_data=True,
                lumi_csv=lumicsv,
                lumi_json=lumijson,
            ))
        else:
            narf_info.update(dict(
                xsec=info["xsec"],
                group=info["group"] if "group" in info else None,
                )
            )
        narf_datasets.append(narf.Dataset(**narf_info))

    # FIXME: the user should probably not be allowed to pass a filter, it is too generic and it requires the user to know
    #        the objects and their attributes on which the filter would act, but the user should not need to know about them
    # It is kept for backward compatibility for now (it is used in the histmakers)
    if filt:
        if isinstance(filt, list):
            # TODO: should probably agree on some convention (e.g. always define a valid group name,
            #       which could become same as name when not specified)
            # next line should only work with groups, but some processes have group=None, so in that case the name is used,
            # and when the name is used the filter has to match (a subset of) the name
            # FIXME: this is a problem if I want to filter data as "Data" (which is the group name), since the sample name is dataPostVFP and the group is None, so it won't match
            # solution 1: always define a group name:
            # solution 2: use the current option passing also the original data set name 
            narf_datasets = list(filter(lambda x: x.group in filt if x.group is not None else any(f in x.name for f in filt), narf_datasets))
        else:
            narf_datasets = list(filter(filt, narf_datasets))
    if excludeGroup:
        if isinstance(excludeGroup, list):
            # some datasets dictionary might not have the group key, hence the narf dataset is defined with group=None
            narf_datasets = list(filter(lambda x: x.group not in excludeGroup if x.group is not None else 1, narf_datasets))
        else:
            narf_datasets = list(filter(excludeGroup, narf_datasets))
            
    for sample in narf_datasets:
        if not sample.filepaths:
            logger.warning(f"Failed to find any files for sample {sample.name}!")

    return narf_datasets

def buildXrdFileList(path, xrd):
    xrdpath = path[path.find('/store'):]
    logger.debug(f"Looking for path {xrdpath}")
    # xrdfs doesn't like wildcards, just use the mount if they are included
    if "*" not in path:
        f = subprocess.check_output(['xrdfs', f'root://{xrd}', 'ls', xrdpath]).decode(sys.stdout.encoding)
        return filter(lambda x: "root" in x[-4:], f.split())
    else:
        return [f"root://{xrd}/{f}" for f in glob.glob(path)]

