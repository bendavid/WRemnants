import ROOT
import pandas as pd
import numpy as np
from scipy.stats import chi2

from utilities import logging, common
from utilities.io_tools import output_tools, tex_tools

import pdb


translate = {
    "asimov": "Asimov",
    "data": "Data",
    "uncorr": "MiNNLO",
    "dyturbo": "DYTurbo",
    "dyturboN3LLp": "DYTurbo N$^{3}$LL' (1D)",
    "dyturboN3LLp2d": "DYTurbo N$^{3}$LL' (2D)",
    "matrix_radish": "MATRIX$+$RadISH",
    "scetlibNP": "SCETLib NP",
    "scetlibN4LL": "SCETLib N$^{4}$LL"
}

def get_chi2(filename):
    try:
        rfile = ROOT.TFile.Open(filename)
        ttree = rfile.Get("fitresults")
        ttree.GetEntry(0)
        val = 2*(ttree.nllvalfull - ttree.satnllvalfull)
        ndf = rfile.Get("obs;1").GetNbinsX() - ttree.ndofpartial
        p = (1-chi2.cdf(val, ndf))*100

        status = ttree.status
        errstatus = ttree.errstatus
        edmval = ttree.edmval

    except IOError as e:
        return 0, 1, 0, 0, 0, 0
        
    return val, ndf, p, status, errstatus, edmval

parser = common.plot_parser()
parser.add_argument("inputs", nargs="+", type=str, help="Paths to fitresult files")
args = parser.parse_args()

logger = logging.setup_logger(__file__, args.verbose, args.noColorLogger)
outdir = output_tools.make_plot_dir(args.outpath, args.outfolder)

df = pd.DataFrame(args.inputs, columns=["path"])

df[["chi2" ,"ndf", "pvalue", "status", "errstatus", "edmval"]] = df["path"].apply(get_chi2).apply(pd.Series)

df["name_parts"] = df["path"].apply(lambda x: [y for y in filter(lambda z: z, x.split("/"))])
df["axes"] = df["name_parts"].apply(lambda x: x[-2].split("_")[1:])
df["channel"] = df["name_parts"].apply(lambda x: x[-2].split("_")[0])
df["column_name"] = df["axes"].apply(lambda x: "-".join(x)) 
df["column_name_ndf"] = df["column_name"] + df["ndf"].apply(lambda x: f" ({x})")

df["dataset"] = df["name_parts"].apply(lambda x: "_".join(x[-1].split("_")[2:]).replace(".root",""))

df["dataset"] = df["dataset"].apply(lambda x: translate.get(x.replace("Corr",""), x.replace("Corr","")))

caption = """
Resulting $\chi^2$ values (and p-values) from fits on different data, and pseudodata sets for the dileptonic analysis.
"""

for channel, df_c in df.groupby("channel"):
    tex_tools.make_latex_table(df_c, output_dir=outdir, output_name=f"table_{channel}", 
        column_title="Axes (bins)", caption=caption, label="Pseudodata", sublabel="",
        column_name="column_name_ndf", row_name="dataset", 
        cell_columns=["chi2", "pvalue"], color_condition=lambda x, y: y < 5, cell_format=lambda x, y: f"${round(x,1)} ({round(y,1)}\%)$")

    tex_tools.make_latex_table(df_c, output_dir=outdir, output_name=f"table_status_{channel}", 
        column_title="Axes", caption="Fit status and error status.", label="Pseudodata", sublabel="",
        column_name="column_name", row_name="dataset", 
        cell_columns=["status", "errstatus"], color_condition=lambda x, y: y !=0, cell_format=lambda x, y: f"{int(x)} ({int(y)})")

    tex_tools.make_latex_table(df_c, output_dir=outdir, output_name=f"table_edmval_{channel}", 
        column_title="Axes", caption="Estimated distance to minimum.", label="Pseudodata", sublabel="",
        column_name="column_name", row_name="dataset", 
        cell_columns=["edmval"], color_condition=lambda x: False, cell_format=lambda x: x)

