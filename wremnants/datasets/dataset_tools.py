import narf
from utilities import logging
import subprocess
import sys
import os
import glob
import random
import pathlib
import socket
#set the debug level for logging incase of full printout 
from wremnants.datasets.datasetDict_v9 import dataDictV9
from wremnants.datasets.datasetDict_v8 import dataDictV8
from wremnants.datasets.datasetDict_gen import genDataDict
from wremnants.datasets.datasetDict_lowPU import dataDictLowPU
import ROOT

logger = logging.child_logger(__name__)

default_nfiles = {
    'WminusmunuPostVFP' : 1700,
    'WplusmunuPostVFP' : 2000,
    'WminustaunuPostVFP' : 400,
    'WplustaunuPostVFP' : 500,
    'ZmumuPostVFP' : 900,
    'ZtautauPostVFP' : 1200,
}

def buildFileListPosix(path):
    # roots, dirs, fnamess = os.walk(path)

    outfiles = []
    for root, dirs, fnames in os.walk(path):
        for fname in fnames:
            if fname.endswith(".root"):
                # print(root, fname)
                outfiles.append(f"{root}/{fname}")

    # print(path, outfiles)

    return outfiles

def buildFileListXrd(path, num_clients = 16):
    xrdprefix = "root://"
    if not path.startswith(xrdprefix):
        raise ValueError(f"Invalid xrootd path {path}")

    path = path[len(xrdprefix):]
    sepidx = path.find("/")
    xrd = path[:sepidx]
    xrdpath = path[sepidx:]


    res = subprocess.run(["xrdfs", f"{xrdprefix}{xrd}", "ls", "-R", xrdpath], capture_output=True)
    xrdfiles = res.stdout.decode(sys.stdout.encoding).split()

    outfiles = []
    for xrdfile in xrdfiles:
        if not xrdfile.endswith(".root"):
            continue

        # construct client string if necessary to force multiple xrootd connections
        # (needed for good performance when a single or small number of xrootd servers is used)
        if num_clients > 1:
            client = f"user_{random.randrange(num_clients)}@"
        else:
            client = ""

        outfiles.append(f"{xrdprefix}{client}{xrd}/{xrdfile}")

    return outfiles

def buildFileList(path):
    xrdprefix = "root://"
    return buildFileListXrd(path) if path.startswith(xrdprefix) else buildFileListPosix(path)

#TODO add the rest of the samples!
def makeFilelist(paths, maxFiles=-1, base_path=None, nano_prod_tags=None, is_data=False, oneMCfileEveryN=None):
    filelist = []
    for orig_path in paths:
        if maxFiles > 0 and len(filelist) >= maxFiles:
            break
        # try each tag in order until files are found
        fallback = False
        for prod_tag in nano_prod_tags:
            format_args=dict(BASE_PATH=base_path, NANO_PROD_TAG=prod_tag)

            path = orig_path.format(**format_args)
            logger.debug(f"Reading files from path {path}")

            files = buildFileList(path)

            if len(files) == 0:
                fallback = True
                logger.warning(f"Did not find any files for tag {prod_tag} matching path {path}!")
            else:
                if fallback:
                    logger.warning(f"Falling back to tag {prod_tag} with path {path}")
                break

        filelist.extend(files)

    if oneMCfileEveryN != None and not is_data:
        tmplist = []
        for i,f in enumerate(filelist):
            if i % oneMCfileEveryN == 0:
                tmplist.append(f)
        logger.warning(f"Using {len(tmplist)} files instead of {len(filelist)}")
        filelist = tmplist

    return filelist if maxFiles < 0 or len(filelist) < maxFiles else random.Random(1).sample(filelist, maxFiles)

def selectProc(selection, datasets):
    if any(selection == x.group for x in datasets):
        # if the selection matches any of the group names in the given dataset, the selection is applied to groups
        return list(filter(lambda x, s=selection: x.group is not None and x.group == s, datasets))
    else:
        # otherwise, the selection is applied to sample names
        return list(filter(lambda x, s=selection: s in x.name, datasets))

def selectProcs(selections, datasets):
    new_datasets = []
    for selection in selections:
        new_datasets += selectProc(selection, datasets)

    # remove duplicates selected by multiple filters
    new_datasets = list(set(new_datasets))
    return new_datasets

def filterProcs(filters, datasets):
    if filters:
        if isinstance(filters, list):
            new_datasets = selectProcs(filters, datasets)
        elif isinstance(filters, str):
            new_datasets = selectProc(filters, datasets)
        else:
            new_datasets = list(filter(filters, datasets))
    else:
        return datasets

    if len(new_datasets) == 0:
        logger.warning("Try to filter processes/groups but didn't find any match. Continue without filtering.")
        return datasets

    return new_datasets

def excludeProcs(excludes, datasets):
    if excludes:
        if isinstance(excludes, list):
            # remove selected datasets
            return list(filter(lambda x: x not in selectProcs(excludes, datasets), datasets))
        elif isinstance(excludes, str):
            # remove selected datasets
            return list(filter(lambda x: x not in selectProc(excludes, datasets), datasets))
        else:
            return list(filter(excludes, datasets))
    else:
        return datasets

def getDataPath(mode=None):
    import socket
    hostname = socket.gethostname()

    if hostname == "lxplus8s10.cern.ch":
        base_path = "/scratch/shared/NanoAOD"
    if hostname == "cmswmass2.cern.ch":
        base_path = "/data/shared/NanoAOD"
    elif "mit.edu" in hostname:
        base_path = "root://submit30.mit.edu//cms/wmass/NanoAOD"
        # base_path = "/scratch/submit/cms/wmass/NanoAOD"
    elif hostname == "cmsanalysis.pi.infn.it":
        base_path = "/scratchnvme/wmass/NANOV9/postVFP"

    # NOTE: If anyone wants to run this at Pisa they'd probably want a different path
    if mode and "lowpu" in mode:
        base_path = f"{base_path}/LowPU/"

    return base_path

def is_zombie(file_path):
    # Try opening the ROOT file and check if it's a zombie file
    file = ROOT.TFile.Open(file_path)
    if not file or file.IsZombie():
        logger.warning(f"Found zombie file: {file_path}")
        return True
    file.Close()
    return False

def getDatasets(maxFiles=default_nfiles, filt=None, excl=None, mode=None, base_path=None, nanoVersion="v9",
                data_tags=["TrackFitV722_NanoProdv2"],
                mc_tags=["TrackFitV722_NanoProdv3", "TrackFitV722_NanoProdv2", "TrackFitV718_NanoProdv1", "0"], oneMCfileEveryN=None, checkFileForZombie=False):
    if maxFiles is None:
        maxFiles=default_nfiles

    if not base_path:
        base_path = getDataPath(mode)
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
    elif mode and "lowpu" in mode:
        dataDict = dataDictLowPU

    narf_datasets = []
    for sample,info in dataDict.items():
        if sample in genDataDict:
            base_path = base_path.replace("NanoAOD", "NanoGen")

        is_data = info.get("group","") == "Data"

        prod_tags = data_tags if is_data else mc_tags
        nfiles = maxFiles
        if type(maxFiles) == dict:
            nfiles = maxFiles[sample] if sample in maxFiles else -1
        paths = makeFilelist(info["filepaths"], nfiles, base_path=base_path, nano_prod_tags=prod_tags, is_data=is_data, oneMCfileEveryN=oneMCfileEveryN)
            
        if checkFileForZombie:
            paths = [p for p in paths if not is_zombie(p)]

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
                lumi_csv=info["lumicsv"],
                lumi_json=info["lumijson"],
                group=info["group"] if "group" in info else None,
            ))
        else:
            narf_info.update(dict(
                xsec=info["xsec"],
                group=info["group"] if "group" in info else None,
                )
            )
        narf_datasets.append(narf.Dataset(**narf_info))
    
    narf_datasets = filterProcs(filt, narf_datasets)
    narf_datasets = excludeProcs(excl, narf_datasets)

    for sample in narf_datasets:
        if not sample.filepaths:
            logger.warning(f"Failed to find any files for sample {sample.name}!")

    return narf_datasets
