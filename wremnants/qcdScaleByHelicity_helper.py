import uproot
import ROOT
import pathlib
import hist
import narf
from .correctionsTensor_helper import makeCorrectionsTensor
from .theory_tools import scale_tensor_axes

data_dir = f"{pathlib.Path(__file__).parent}/data/"

def makeQCDScaleByHelicityHelper(input_path=f"{data_dir}/angularCoefficients"):
    axis_chargeVgen = hist.axis.Regular(
        2, -2, 2, name="chargeVgen", underflow=False, overflow=False
    )
    axis_helicity = hist.axis.Integer(
        -1, 9, name="helicity", overflow=False, underflow=False
    )
    f = uproot.open(f"{input_path}/fractions_minus_2022.root")
    nom = f["unpol_minus_nominal"].to_hist()
    axis_yVgen, axis_ptVgen = [hist.axis.Variable(ax.edges, name=name) for ax,name in zip(nom.axes, ["yVgen", "ptVgen"])]
    corrh = hist.Hist(axis_yVgen, axis_ptVgen, axis_chargeVgen, axis_helicity, *scale_tensor_axes)

    for i,charge in enumerate(["minus", "plus"]):
        f = uproot.open(f"{input_path}/fractions_{charge}_2022.root")
        for k,coeff in enumerate(["unpol"] + [f"a{m}" for m in range(8)]):
            corrh.view(flow=True)[...,i, k, 1, 1] = f[f"{coeff}_{charge}_nominal"].to_hist().values(flow=True)
            for j,muR in enumerate(["muRDown", "", "muRUp"]):
                for l,muF in enumerate(["muFDown", "", "muFUp"]):
                    if muR != "" and muF == "":
                        name = muR
                    elif muF != "" and muR == "":
                        name = muF
                    elif "Up" in muR and "Up" in muF:
                        name = "muRmuFUp"
                    elif "Down" in muR and "Down" in muF:
                        name = "muRmuFDown"
                    else:
                        continue

                    rthist = f[f"{coeff}_{charge}_{name}"].to_hist()
                    corrh.view(flow=True)[...,i, k, j, l] = rthist.values(flow=True)

    #should be uniformly 1.0 for 1+cos^2 theta term
    corrh.values(flow=True)[...,0, :, :] = 1.0
    corrh.variances(flow=True)[...,0, :, :] = 0.0

    return makeCorrectionsTensor(corrh, ROOT.wrem.QCDScaleByHelicityCorrectionsHelper, tensor_rank=3)
