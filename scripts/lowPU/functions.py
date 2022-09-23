
import sys,array,math,os,copy,shutil,decimal


import ROOT
ROOT.gROOT.SetBatch(True)
ROOT.gStyle.SetOptStat(0)
ROOT.gStyle.SetOptTitle(0)


def prepareDir(outDir, remove=True):

    if os.path.exists(outDir) and os.path.isdir(outDir) and remove: shutil.rmtree(outDir)
    os.system("mkdir -p %s" % outDir)
    os.system("cp /eos/user/j/jaeyserm/www/wmass/index.php %s" % outDir)

    
def doOverlow(h):

    n = h.GetNbinsX()
    h.SetBinContent(1, h.GetBinContent(0) + h.GetBinContent(1))
    h.SetBinContent(n, h.GetBinContent(n+1) + h.GetBinContent(n))
    h.SetBinError(1, math.hypot(h.GetBinError(0), h.GetBinError(1)))
    h.SetBinError(n, math.hypot(h.GetBinError(n+1), h.GetBinError(n)))
    h.SetBinContent(0, 0)
    h.SetBinContent(n+1, 0)
    h.SetBinContent(0, 0)
    h.SetBinContent(n+1, 0)
    
    return h    

def parseProc(groups, histCfg, procName, syst="", rebin=1):

    axis = histCfg['axis']
    hNames = histCfg['name'].split(",")

    
    bhist = None
    for hName in hNames:
        
        bhist = groups.readProc(hName, procName, axis=axis)
        if bhist == None: continue
        label = "%s_%s" % (hName, procName)
        break
        
    print(bhist)
    rhist = narf.hist_to_root(bhist)
    rhist.Rebin(rebin)
    rhist.SetName(label)
    rhist = doOverlow(rhist)

    print("Get histogram %s, yield=%d" % (label, rhist.Integral()))
    return rhist
    
    

def Rebin(h, newbins, binWidth=True):

    if isinstance(newbins, int):
        h.Rebin(newbins)
        if binWidth: h.Scale(1, "width")
        return h
    else:
        mybins = array.array('d', newbins)
        h1 = h.Rebin(len(mybins)-1, h.GetName(), mybins)
        if binWidth: h1.Scale(1, "width")
        return h1


def drange(x, y, jump):
    while x < y:
        yield float(x)
        #x += decimal.Decimal(jump)
        x += jump
        