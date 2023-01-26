import string
import os
import sys 
import subprocess
import datetime
import time
import pickle
import lz4.frame
import logging
from utilities import common

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

def metaInfoDict(exclude_diff='notebooks'):
    meta_data = {"time" : str(datetime.datetime.now()), "command" : ' '.join(sys.argv)}
    if subprocess.call(["git", "branch"], stderr=subprocess.STDOUT, stdout=open(os.devnull, 'w')) != 0:
        meta_data["git_info"] = {"hash" : "Not a git repository!",
                "diff" : "Not a git repository"}
    else:
        meta_data["git_hash"] = subprocess.check_output(['git', 'log', '-1', '--format="%H"'], encoding='UTF-8')
        diff_comm = ['git', 'diff']
        if exclude_diff:
            diff_comm.extend(['--', f":!{exclude_diff}"])
        meta_data["git_diff"] = subprocess.check_output(diff_comm, encoding='UTF-8')

    return meta_data

def analysis_debug_output(results):
    logging.debug("")
    logging.debug("Unweighted events (before cut)")
    logging.debug("-"*30)
    for key,val in results.items():
        if "event_count" in val:
            logging.debug(f"Dataset {key.ljust(30)}:  {val['event_count']}")
            logging.debug("-"*30)
    logging.debug("")

def writeMetaInfoToRootFile(rtfile, exclude_diff='notebooks'):
    import ROOT
    meta_dict = metaInfoDict(exclude_diff)
    d = rtfile.mkdir("meta_info")
    d.cd()
    
    for key, value in meta_dict.items():
        out = ROOT.TNamed(str(key), str(value))
        out.Write()

def write_analysis_output(results, outfile, args):
    analysis_debug_output(results)
    results.update({"meta_info" : metaInfoDict()})

    to_append = []
    if args.theory_corr and not args.theory_corr_alt_only:
        to_append.append(args.theory_corr[0]+"Corr")
    if hasattr(args, "uncertainty_hist") and args.uncertainty_hist != "nominal":
        to_append.append(args.uncertainty_hist)
    if args.postfix:
        to_append.append(args.postfix)
    if args.maxFiles > 0:
        to_append.append(f"maxFiles{args.maxFiles}")

    if to_append:
        outfile = outfile.replace(".pkl.lz4", f"_{'_'.join(to_append)}.pkl.lz4")

    output = os.path.join(args.outfolder, outfile)

    time0 = time.time()
    print(f"writing output file {output} ...")
    with lz4.frame.open(output, "wb") as f:
        pickle.dump(results, f, protocol = pickle.HIGHEST_PROTOCOL)
    print("writing output:", time.time()-time0)
