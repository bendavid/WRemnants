import uproot
import hist
import argparse
from wremnants import theory_corrections,plot_tools,theory_tools
from wremnants import boostHistHelpers as hh
import pickle
import lz4.frame
import os
from matplotlib import cm
import matplotlib.pyplot as plt
import logging

logging.basicConfig(level=logging.INFO)

s = hist.tag.Slicer()
# Map hist_name argument to the actual hist and it's axis stored in files
lookup = {
    "minnlo" : {
        "dirs" : {
            "z" : ("ZmumuPostVFP", "output"),
            "wp" : ("WplusmunuPostVFP", "output"),
            "wm" : ("WminusmunuPostVFP", "output"),
        },
        "ptV" : {
            "hist" : "nominal_gen",
            "axis" : "ptVgen",
            "action" : lambda x: x[{"massVgen" : s[0:x.axes["massVgen"].size:hist.sum]}],
        },
        "absYV" : {
            "hist" : "nominal_gen",
            "axis" : "absYVgen",
            "action" : lambda x: x[{"massVgen" : s[0:x.axes["massVgen"].size:hist.sum], "ptVgen" : s[0:40.j:hist.sum]}],
        },
        "sigma4_ptV" : {
            "hist" : "helicity_moments_scale",
            "axis" : "ptVgen",
            "action" : lambda x: theory_tools.scale_angular_moments(
                x[{"muRfact" : 1.j, "muFfact" : 1.j, "massVgen" : s[0:x.axes["massVgen"].size:hist.sum]}])[{"helicity" : 4.j}],
        },
        "sigma4_absYV" : {
            "hist" : "helicity_moments_scale",
            "axis" : "absYVgen",
            "action" : lambda x: theory_tools.scale_angular_moments(
                x[{"muRfact" : 1.j, "muFfact" : 1.j, "massVgen" : s[0:x.axes["massVgen"].size:hist.sum], "ptVgen" : s[0:40.j:hist.sum]}])[{"helicity" : 4.j}],
        }
    },
    "scetlib" : {
        "dirs" : {"z" : ("W"), "wp" : ("W"), "wm" : ("W")},
        "ptV" : {
            "hist" : "inclusive_pT_y_m_scetlib",
            "axis" : "pt",
        },
        "absYV" : {
            "hist" : "inclusive_pT_y_m_scetlib",
            "axis" : "absy",
            "action" : lambda x: hh.makeAbsHist(x[{"pt" : s[0:40.j:hist.sum]}], "y"),
        },
        "sigma4_ptV" : {
            "axis" : "pt",
        },
        "sigma4_absYV" : {
            "axis" : "absy",
            "action" : lambda x: hh.makeAbsHist(x[{"pt" : s[0:40.j:hist.sum]}], "y"),
        },
    },
    "matrix_radish" : {
        "ptV" : {
            "z": {
                "axis" : "xaxis",
            },
            "wm" : {
                "axis" : "xaxis",
            },
        },
        #"absYV" : {
        #    "hist" : 
    },
    "dyturbo" : {
        "ptV" : {
            "axis" : "ptV",
        }
    },
}

parser = argparse.ArgumentParser()
parser.add_argument("-s", "--scetlib_files", default=[], type=str, nargs='+', help="SCETlib input files (.npz format)")
parser.add_argument("-m", "--minnlo_files", default=[], type=str, nargs='+', help="minnlo input files (.pkl format)")
parser.add_argument("-mr", "--matrix_radish_files", default=[], type=str, nargs='+', help="MATRIX+RadISH input files (.dat format)")
parser.add_argument("-dt", "--dyturbo_files", default=[], type=str, nargs='+', help="DYTurbo input files (.txt format)")
parser.add_argument("-n", "--hist_name", type=str, choices=lookup["minnlo"].keys(), help="Observable to plot", required=True)
parser.add_argument("-p", "--proc", type=str, choices=lookup["minnlo"]["dirs"].keys(), help="Process", required=True)
parser.add_argument("-r", "--rrange", type=float, nargs=2, default=[0.9, 1.1], help="y range for ratio plot")
parser.add_argument("--ylim", type=float, nargs=2, help="Range for y axis (if not specified, range set automatically)")
parser.add_argument("-o", "--outpath", type=str, default=os.path.expanduser("~/www/WMassAnalysis"), help="Base path for output")
parser.add_argument("-f", "--outfolder", type=str, default="test", help="Subfolder for output")
parser.add_argument("-a", "--name_append", type=str, help="Name to append to file name")
parser.add_argument("--ratio_ref", type=str, default="minnlo", choices=["minnlo", "scetlib", "matrix_radish"], help="Reference hist for ratio")
parser.add_argument("--keep_full_range", action='store_true', help="Store the full range of all hists, even if it exceeds other hist ranges")
args = parser.parse_args()

cmap = cm.get_cmap("tab10")
lookup["minnlo"]["colors"] = ["red"]+[cmap(i) for i in range(len(args.minnlo_files)-1)]
lookup["scetlib"]["colors"] = ["purple"]+[cmap(9-i) for i in range(len(args.scetlib_files)-1)]
lookup["matrix_radish"]["colors"] = ["green"]+[cmap(i+len(args.minnlo_files)) for i in range(len(args.scetlib_files)-1)]
lookup["dyturbo"]["colors"] = ["blue"]+[cmap(i) for i in range(len(args.dyturbo_files))]

xlabels = {
    "ptV" : r"p$_{T}^{%s}$ (GeV)" % ("W" if "w" in args.proc else "Z"), 
    "absYV" : r"$|\mathrm{y}^{%s}|$"  % ("W" if "w" in args.proc else "Z"),
}
xlabels["sigma4_ptV"] = xlabels["ptV"]
xlabels["sigma4_absYV"] = xlabels["absYV"]

ylabels = {
    "sigma4_ptV" : "$\sigma_{4}$/bin",
    "sigma4_absYV" : "$\sigma_{4}$/bin",
    "a4_ptV" : "$\A_{4}$/bin",
    "a4_absYV" : "$\A_{4}$/bin",
}

def transform_and_project(histND, scale, axis_name, action=None):
    if action is not None:
        histND = action(histND)
    hist1D = histND.project(axis_name)*scale
    return hist1D

def read_pickle_hist(proc, file_name, lookup, hist_name):
    with lz4.frame.open(file_name) as f:
        output = pickle.load(f)

    #TODO: This code is so terrible, fix it!
    levels = {proc : []} if "dirs" not in lookup else lookup["dirs"]
    hist_out = output
    scale = 1.
    for level in levels[proc]:
        if "dataset" in hist_out:
            scale *= hist_out["dataset"]["xsec"]
        if "weight_sum" in hist_out:
            scale = scale/hist_out["weight_sum"]
        hist_out = hist_out[level]

    if hist_name not in lookup:
        raise ValueError(f"Couldn't find hist {hist_name} in file {file_name}")

    hist_info = lookup[hist_name]
    histND = hist_out[hist_info["hist"]]
    action = hist_info["action"] if "action" in hist_info else None
    return transform_and_project(histND, scale, hist_info["axis"], action)

def read_scetlib_hist(proc, file_name, lookup, hist_name):
    charge = 0 if proc == "z" else (-1 if proc == "wp" else 1)
    nonsing = file_name.replace(".npz", "_nons.npz")
    if not os.path.isfile(nonsing):
        logging.warning("Didn't find the non-singular contribution. Will make comparisons without it")
        nonsing = ""
    histND = theory_corrections.read_scetlib_hist(file_name, nonsing=nonsing, charge=charge, flip_y_sign="A4" in file_name)

    hist_info = lookup[hist_name]
    action = hist_info["action"] if "action" in hist_info else None
    return transform_and_project(histND, 1.0, hist_info["axis"], action)

def read_matrix_radish_hist(proc, filename, lookup, hist_name):
    hist_info = lookup[hist_name][proc]
    histND = theory_corrections.read_matrixRadish_hist(filename, hist_info["axis"])
    scale = 1.0 if "scale" not in hist_info else hist_info["scale"]

    action = hist_info["action"] if "action" in hist_info else None
    return transform_and_project(histND, scale, hist_info["axis"], action)

def read_hists(proc, files, lookup, hist_name):
    hists = []
    for fname in files:
        if ".pkl.lz4" in fname[-8:]:
            func = read_pickle_hist
        if ".npz" in fname[-4:]:
            func = read_scetlib_hist
        if ".dat" in fname[-4:]:
            func = read_matrix_radish_hist
        if ".txt" in fname[-4:]:
            func = read_dyturbo_hist
        hists.append(func(proc, fname, lookup, hist_name))
    return hists

if not args.scetlib_files and not args.minnlo_files:
    raise ValueError("Must specify at least one filename")

outdir = plot_tools.make_plot_dir(args.outpath, args.outfolder)

all_hists = []
all_labels = []
all_colors = []

generators_info = [
	("minnlo", "MiNNLO (NNLO+PS)"), 
	("scetlib", "SCETlib (N$^{3}$LL)"),
	("matrix_radish", "MATRIX+RadISH (NNLO+N$^{3}$LL)"),
	("dyturbo", "DYTurbo (NNLO+N$^{3}$LL)"),
]

short_name = {
    "minnlo" : "MiNNLO",
    "scetlib" : "SCETlib",
    "matrix_radish" : "MAT+Rad",
    "dyturbo" : "DYTurbo",
}


generators_info.insert(0, generators_info.pop([x[0] for x in generators_info].index(args.ratio_ref)))

for generator, label in generators_info:
    files = getattr(args, f"{generator}_files")
    if files:
        info = lookup[generator]
        if args.hist_name not in info:
            continue
        if generator == "dyturbo":
            hists = [theory_corrections.read_dyturbo_hist(files, axis=info[args.hist_name]["axis"])]
        else:
            hists = read_hists(args.proc, files, info, args.hist_name)
        all_hists.extend(hists)
        all_colors.extend([info["colors"][c] for c in range(len(hists))])
        all_labels.append(label)
        if len(hists) > 1:
            all_labels.extend([f"{label} (alt {i})" for i in range(1, len(hists[1:]+1))])

if len(all_hists) > 1:
    all_hists = hh.rebinHistsToCommon(all_hists, axis_idx=0, keep_full_range=args.keep_full_range)

ylabel = "$\sigma$/bin" if args.hist_name not in ylabels else ylabels[args.hist_name]
fig = plot_tools.makePlotWithRatioToRef(all_hists, colors=all_colors, labels=all_labels, alpha=0.7, ylim=args.ylim,
        rrange=args.rrange, ylabel=ylabel, xlabel=xlabels[args.hist_name], rlabel=f"x/{short_name[args.ratio_ref]}", binwnorm=1.0, nlegcols=1)

outname = f"TheoryCompHist_{args.proc}_{args.hist_name}" + ("_"+args.name_append if args.name_append else "")
plot_tools.save_pdf_and_png(outdir, outname)
plot_tools.write_index_and_log(outdir, outname)
