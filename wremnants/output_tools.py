import string
import os
import sys 
import subprocess
import datetime

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

def metaInfoDict():
    meta_data = {"time" : str(datetime.datetime.now()), "command" : ' '.join(sys.argv)}
    if subprocess.call(["git", "branch"], stderr=subprocess.STDOUT, stdout=open(os.devnull, 'w')) != 0:
        meta_data["git_info"] = {"hash" : "Not a git repository!",
                "diff" : "Not a git repository"}
    else:
        meta_data["git_hash"] = subprocess.check_output(['git', 'log', '-1', '--format="%H"'], encoding='UTF-8')
        meta_data["git_diff"] = subprocess.check_output(['git', 'diff'], encoding='UTF-8')

    return meta_data
