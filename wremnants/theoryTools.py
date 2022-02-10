def qcdByHelicityLabels():
    coeffs = ["const"]+[f"a{i}" for i in range(8)]
    scaleVars = ["muRmuF", "muR", "muF"]
    return  [f"{var}_{coeff}{t}" for var in scaleVars for t in ["Up", "Down"] for coeff in coeffs]

def qcdScaleNames():
    # Exclude central and extreme variations
    shifts = ["muRmuFDown", "muRDown", "", "muFDown", "", "muFUp", "muRUp", "", "muRmuFUp"]
    return ["_".join(["QCDscale", s]) if s != "" else s for s in shifts]

def massWeightNames(matches=None):
    central=10
    nweights=21
    names = [f"massShift{int(abs(central-i)*10)}MeV{'Down' if i < central else 'Up'}" for i in range(nweights)]
    # If name is "" it won't be stored
    return [x if not matches or any(y in x for y in matches) else "" for x in names]

def pdfNames(cardTool, pdf, skipFirst=False):
    size = 103
    names = cardTool.mirrorNames(f"pdf{{i}}{pdf}", size)
    if skipFirst:
        names[0] = ""
        names[size] = ""
    if pdf == "NNPDF31":
        names[size-2] = "pdfAlphas002Up"
        names[size-1] = "pdfAlphas002Down"
        # Drop the mirrored alphaS variations
        names[size*2-2] = ""
        names[size*2-1] = ""
    return names


