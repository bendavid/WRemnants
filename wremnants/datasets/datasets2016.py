import narf
import logging
import subprocess
import glob
import pathlib
import socket
import logging
#set the debug level for logging incase of full printout 
from wremnants.datasets.datasetDict_v9 import dataDictV9, dataDictV9_pisa
from wremnants.datasets.datasetDict_v8 import *

lumicsv = f"{pathlib.Path(__file__).parent.parent}/data/bylsoutput.csv"
lumijson = f"{pathlib.Path(__file__).parent.parent}/data/Cert_271036-284044_13TeV_Legacy2016_Collisions16_JSON.txt"

#TODO add the rest of the samples!
def makeFilelist(paths, maxFiles=-1):
    filelist = []
    for path in paths:
        filelist.extend(glob.glob(path) if path[:4] != "/eos" else buildXrdFileList(path, "eoscms.cern.ch"))
    return filelist if maxFiles < 0 else filelist[:maxFiles]

def getNarfDataset(sampleName, maxFiles, sampleDict, isData, isWorZ=True):
    if isData:
        logging.debug('Sample ', sampleName, ' read from : ', sampleDict[sampleName]["filepaths"])
        nData = narf.Dataset(name = sampleDict[sampleName]["name"],
                                   filepaths = makeFilelist(sampleDict[sampleName]["filepaths"], maxFiles),
                                   is_data = True,
                                   lumi_csv = lumicsv,
                                   lumi_json = lumijson
        )
        return nData
    else:
        logging.debug('Sample ', sampleName, ' read from : ', sampleDict[sampleName]["filepaths"])
        nMC = narf.Dataset(name = sampleDict[sampleName]["name"],
                           filepaths = makeFilelist(sampleDict[sampleName]["filepaths"], maxFiles),
                           is_data = False,
                           xsec = sampleDict[sampleName]['xsec'],
        )
        if not isWorZ: nMC.group = sampleDict[sampleName]['group']
        return nMC
        
def getDatasets(maxFiles=-1, filt=None, mode=None, nanoVersion = "v9"):
    dataDict = dataDictV9_pisa if socket.gethostname() == 'cmsanalysis.pi.infn.it' else dataDictV9
    if nanoVersion != "v9":
        dataDict = dataDictV8
        print('Using data dict V8')

    dataPostVFP = getNarfDataset("dataPostVFP", maxFiles, dataDict, True)

    ZmmPostVFP = getNarfDataset("ZmmPostVFP", maxFiles, dataDict, False, True)

    ZttPostVFP = getNarfDataset("ZttPostVFP", maxFiles, dataDict, False, True)

    WpmunuPostVFP = getNarfDataset("WpmunuPostVFP", maxFiles, dataDict, False, True)

    WmmunuPostVFP = getNarfDataset("WmmunuPostVFP", maxFiles, dataDict, False, True)

    WptaunuPostVFP = getNarfDataset("WptaunuPostVFP", maxFiles, dataDict, False, True)

    WmtaunuPostVFP = getNarfDataset("WmtaunuPostVFP", maxFiles, dataDict, False, True)

    ttbarlnuPostVFP = getNarfDataset("ttbarlnuPostVFP", maxFiles, dataDict, False, False)

    ttbarlqPostVFP = getNarfDataset("ttbarlqPostVFP", maxFiles, dataDict, False, False)

    singleTop_schanLepDecaysPostVFP = getNarfDataset("singleTop_schanLepDecaysPostVFP", maxFiles, dataDict, False, False)

    singleTop_tWAntitopPostVFP = getNarfDataset("singleTop_tWAntitopPostVFP", maxFiles, dataDict, False, False)

    singleTop_tchanAntitopPostVFP = getNarfDataset("singleTop_tchanAntitopPostVFP", maxFiles, dataDict, False, False)

    singleTop_tchanTopPostVFP = getNarfDataset("singleTop_tchanTopPostVFP", maxFiles, dataDict, False, False)

    wwPostVFP = getNarfDataset("wwPostVFP", maxFiles, dataDict, False, False)

    wzPostVFP = getNarfDataset("wzPostVFP", maxFiles, dataDict, False, False)

    zz2l2nuPostVFP = getNarfDataset("zz2l2nuPostVFP", maxFiles, dataDict, False, False)

    allPostVFP = [dataPostVFP,
                  WpmunuPostVFP, WmmunuPostVFP, WptaunuPostVFP, WmtaunuPostVFP,
                  ZmmPostVFP, ZttPostVFP,
                  ttbarlnuPostVFP, ttbarlqPostVFP,
                  singleTop_schanLepDecaysPostVFP, singleTop_tWAntitopPostVFP, singleTop_tchanAntitopPostVFP, singleTop_tchanTopPostVFP,
                  wwPostVFP, wzPostVFP, zz2l2nuPostVFP]
    # ,WmmunuPostVFP_LZ4_4,WmmunuPostVFP_LZMA_9]

    allPostVFP_gen = allPostVFP[1:]

    samples = allPostVFP if mode != "gen" else allPostVFP_gen
    if filt:
        return list(filter(filt, samples))
    else:
        return samples

def buildXrdFileList(path, xrd):
    xrdpath = path[path.find('/store'):]
    logging.debug(f"Looking for path {xrdpath}")
    f = subprocess.check_output(['xrdfs', f'root://{xrd}', 'ls', xrdpath]).decode(sys.stdout.encoding)
    return filter(lambda x: "root" in x[-4:], f.split())
