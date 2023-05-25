import string
import os
import sys 
import subprocess
import datetime
import time
import hdf5plugin
import h5py
import narf
import numpy as np
import re
from utilities import logging

logger = logging.child_logger(__name__)

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

def script_command_to_str(argv, parser_args):
    call_args = np.array(argv[1:], dtype=object)
    match_expr = "|".join(["^-+([a-z]+[1-9]*-*)+"]+([] if not parser_args else [f"^-*{x.replace('_', '.')}" for x in vars(parser_args).keys()]))
    if call_args.size != 0:
        flags = np.vectorize(lambda x: bool(re.match(match_expr, x)))(call_args)
        special_chars = np.vectorize(lambda x: not x.isalnum())(call_args)
        select = ~flags & special_chars
        if np.count_nonzero(select):
            call_args[select] = np.vectorize(lambda x: f"'{x}'")(call_args[select])
    return " ".join([argv[0], *call_args])

def metaInfoDict(exclude_diff='notebooks', args=None):
    meta_data = {
        "time" : str(datetime.datetime.now()), 
        "command" : script_command_to_str(sys.argv, args),
        "args": {a: getattr(args,a) for a in vars(args)}
    }
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
    logger.debug("")
    logger.debug("Unweighted events (before cut)")
    logger.debug("-"*30)
    for key,val in results.items():
        if "event_count" in val:
            logger.debug(f"Dataset {key.ljust(30)}:  {val['event_count']}")
            logger.debug("-"*30)
    logger.debug("")

def writeMetaInfoToRootFile(rtfile, exclude_diff='notebooks', args=None):
    import ROOT
    meta_dict = metaInfoDict(exclude_diff, args=args)
    d = rtfile.mkdir("meta_info")
    d.cd()
    
    for key, value in meta_dict.items():
        out = ROOT.TNamed(str(key), str(value))
        out.Write()

def write_analysis_output(results, outfile, args, update_name=True):
    analysis_debug_output(results)
    results.update({"meta_info" : metaInfoDict(args=args)})

    to_append = []
    if args.theoryCorr and not args.theoryCorrAltOnly:
        to_append.append(args.theoryCorr[0]+"Corr")
    if hasattr(args, "uncertainty_hist") and args.uncertainty_hist != "nominal":
        to_append.append(args.uncertainty_hist)
    if args.postfix:
        to_append.append(args.postfix)
    if args.correlatedNonClosureNP:
        to_append.append("NonClosureCorl")
    if args.maxFiles > 0:
        to_append.append(f"maxFiles{args.maxFiles}")

    if to_append and update_name:
        outfile = outfile.replace(".hdf5", f"_{'_'.join(to_append)}.hdf5")

    if args.outfolder:
        if not os.path.exists(args.outfolder):
            logger.info(f"Creating output folder {args.outfolder}")
            os.makedirs(args.outfolder)
        outfile = os.path.join(args.outfolder, outfile)

    time0 = time.time()
    with h5py.File(outfile, 'w') as f:
        narf.ioutils.pickle_dump_h5py("results", results, f)
    logger.info(f"Writing output: {time.time()-time0}")
    logger.info(f"Output saved in {outfile}")
