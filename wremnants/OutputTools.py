import string
import os

def readTemplate(templateFile, templateDict, filt=None):
    if not os.path.isfile(templateFile):
        raise ValueError("Template file %s is not a valid file!" % templateFile)
    with open(templateFile, "r") as tf:
        lines = filter(filt, tf.readlines()) if filt else tf.readlines()
        source = string.Template("".join(lines))
    filled = source.substitute(templateDict)
    return filled

def fillTemplatedFile(templateFile, outFile, templateDict, append=False):
    filled = readFromTempate(templateFile, templateDict)
    with open(outFile, "w" if not append else "a") as outFile:
        outFile.write(result)


def resetBinError(histo):
    # TODO: move to a C++ helper function to be faster
    if "TH1" in histo.ClassName():
        for b in range(0,histo.GetNbinsX()+2):
            histo.SetBinError(b, 0.0)
    elif "TH2" in histo.ClassName():
        for bx in range(0,histo.GetNbinsX()+2):
            for by in range(0,histo.GetNbinsY()+2):
                histo.SetBinError(bx, by, 0.0)
    elif "TH3" in histo.ClassName():
        for bx in range(0,histo.GetNbinsX()+2):
            for by in range(0,histo.GetNbinsY()+2):
                for bz in range(0,histo.GetNbinsZ()+2):
                    histo.SetBinError(bx, by, bz, 0.0)
