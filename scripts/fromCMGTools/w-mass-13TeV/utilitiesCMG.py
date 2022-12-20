import ROOT, os, sys, re, array, math, json
import numpy as np

class templateBinning:
    def __init__(self,etaBins=[],ptBins=[]):
        self.etaBins = etaBins
        self.ptBins = ptBins
        self.Neta = len(etaBins)-1
        self.Npt  = len(ptBins)-1
        self.NTotBins = self.Neta * self.Npt

    def printBin(self):
        print("###########################")
        print("Binning: eta-pt on x-y axis")
        print("eta bins: %s" % str(self.Neta))
        print("pt  bins: %s" % str(self.Npt))
        print("")

    def printBinAll(self):
        print("###########################")
        print("Binning: eta-pt on x-y axis (%d bins)" % self.NTotBins)
        print("eta bins: %s" % str(self.Neta))
        print("%s" % str(self.etaBins))
        print("-"*20)
        print("pt  bins: %s" % str(self.Npt))
        print("%s" % str(self.ptBins))
        print("-"*20)
        print("")

def getArrayParsingString(inputString, verbose=False, makeFloat=False):
    # convert string [a,b,c,...] to list of a b c ...
    tmp = inputString.replace('[','').replace(']','')
    tmp = tmp.split(',')
    if verbose:
        print("Input:",inputString)
        print("Output:",tmp)
    if makeFloat:
        ret = [float(x) for x in tmp]
    else:
        ret = tmp
    return ret

def getBinning(inputBins, whichBins="reco"):

    # whichBins can be reco or gen
    if whichBins not in ["reco", "gen"]:
        print("Error in function getBinning(): whichBins must be 'reco' or 'gen'. Exit")
        exit()

    # case in which we are passing a file containing the binning and not directly the binning itself
    if inputBins.startswith("file=") or re.match(".*binningPtEta.*.txt",inputBins):
        etaPtbinningFile = inputBins.replace("file=","")
        with open(etaPtbinningFile) as f:
            content = f.readlines()
        for x in content:
            if str(x).startswith(whichBins):
                tmpbinning = (x.split(whichBins+":")[1]).strip()
            else:
                continue
        etabinning = tmpbinning.split('*')[0]    # this is like [a,b,c,...], and is of type string. We nedd to get an array  
        ptbinning  = tmpbinning.split('*')[1]
    else:
        etabinning = inputBins.split('*')[0]    # this is like [a,b,c,...], and is of type string. We nedd to get an array  
        ptbinning  = inputBins.split('*')[1]
    etabinning = list(float(i) for i in etabinning.replace('[','').replace(']','').split(','))
    ptbinning  = list(float(i) for i in ptbinning .replace('[','').replace(']','').split(','))
    binning = [etabinning, ptbinning] 
    return binning


class util:

    def __init__(self):
        #self.cmssw_version = os.environ['CMSSW_VERSION']
        #self.isRecentRelease = (len(self.cmssw_version) and int(self.cmssw_version.split('_')[1]) > 8)
        self.isRecentRelease = True

    def checkHistInFile(self, h, hname, fname, message=""):
        if not h:
            print("Error {msg}: I couldn't find histogram {h} in file {f}".format(msg=message,h=hname,f=fname))
            quit()


    def solvePol2(self, a,b,c):
    
        # calculate the discriminant
        d = (b**2) - (4*a*c)
    
        if not a or d < 0:
            return (0,0,0)
        
        # find two solutions
        sol1 = (-b-math.sqrt(d))/(2*a)
        sol2 = (-b+math.sqrt(d))/(2*a)
    
        bestfit = -1.*b/(2.*a)
    
        return (bestfit, sol1, sol2)
    
    def graphStyle(self, graph, style=20, color=ROOT.kOrange+7, size=1.0, titleY='-2 #Delta ln L', rangeY=(-0.01,4.0) ):
        graph.SetMarkerStyle(style)
        graph.SetMarkerColor(color)
        graph.SetLineWidth  (2)
        graph.SetMarkerSize(size)
        graph.GetYaxis().SetTitle(titleY)
        if rangeY:
            graph.GetYaxis().SetRangeUser(rangeY[0], rangeY[1])
    
    
    def getGraph(self, infile, par, norm, treename='limit'):
        f = ROOT.TFile(infile,'read')
        tree = f.Get(treename)
        vals = []
        normval = norm if norm else 1.
        for ev in tree:
            vals.append( [getattr(ev, par)/normval, 2.*ev.deltaNLL] )
        vals = sorted(vals)
        graph = ROOT.TGraph(len(vals), array.array('d', [x[0] for x in vals]), array.array('d', [y[1] for y in vals]) )
        self.graphStyle(graph)
        graph.GetXaxis().SetTitle(par)
        graph.SetTitle('scan for '+par)
        return graph

    def getErrorFromGraph(self, graph):
        graph.Fit('pol2')
        tmp_fit = graph.GetFunction('pol2')
        (best, sol1, sol2) = self.solvePol2(tmp_fit.GetParameter(2), tmp_fit.GetParameter(1), tmp_fit.GetParameter(0)-1)
        return (best, sol1, sol2)

    def getRebinned(self, ybins, charge, infile, ip):
        histo_file = ROOT.TFile(infile, 'READ')
    
        pstr = 'central' if not ip else 'pdf{ip}'.format(ip=ip)
    
        histos = {}
        for pol in ['left','right','long']:
            cp = '{ch}_{pol}'.format(ch=charge,pol=pol if not pol == 'long' else 'right')
    
            keys = histo_file.GetListOfKeys()
            for k in keys:
                if 'w{ch}'.format(ch=charge) in k.GetName() and pol in k.GetName() and pstr in k.GetName():
                    name = k.GetName()
            histo = histo_file.Get(name)# 'w{ch}_wy_W{ch}_{pol}'.format(ch=charge, pol=pol))
            conts = []
            epsilon = 0.000001 
            # val is a bin boundary, add epsilon to be sure to catch the correct bin (float can be truncated and findBin might return incorrect bin)
            for iv, val in enumerate(ybins[cp][:-1]):
                err = ROOT.Double()
                istart = histo.FindBin(val+epsilon)
                iend   = histo.FindBin(ybins[cp][iv+1]+epsilon)
                val = histo.IntegralAndError(istart, iend-1, err)/36000. ## do not include next bin
                conts.append(float(val))
            histos[pol] = conts
        histo_file.Close()
        return histos

    def getXSecFromShapes(self, ybins, charge, infile, ip, nchannels=1, polarizations = ['left','right','long'], excludeYbins = [], generator='fewz3p1'):
        ## the xsec used in the templates is FEWZ 3.1
        rescale = self.wxsec(generator)[0]/self.wxsec('fewz3p1')[0]
        values = {}
        if not infile:
            for pol in polarizations: 
                cp = '{ch}_{pol}'.format(ch=charge,pol=pol)
                xsecs = []
                for iv in xrange(len(ybins[cp][:-1])):
                    if any(iv == x for x in excludeYbins): continue
                    xsecs.append(0.)
                values[pol] = xsecs
            return values

        histo_file = ROOT.TFile(infile, 'READ')
    
        pstr = '' if not ip else '_pdf{ip}Up'.format(ip=ip)
    
        for pol in polarizations:
            cp = '{ch}_{pol}'.format(ch=charge,pol=pol)
            xsecs = []
            for iv, val in enumerate(ybins[cp][:-1]):                
                if val in excludeYbins: continue
                name = 'x_W{ch}_{pol}_Ybin_{iy}{suffix}'.format(ch=charge,pol=pol,iy=iv,ip=ip,suffix=pstr)                
                histo = histo_file.Get(name)
                val = rescale*histo.Integral()/36000./float(nchannels) # xsec file yields normalized to 36 fb-1
                xsecs.append(float(val))
            values[pol] = xsecs
        histo_file.Close()
        return values

    def getQCDScaleEnvelope(self, ybins, charge, infile, nchannels=1, polarizations = ['left','right','long'], excludeYbins = [], doAlphaS=False, generator='mcatnlo'):
        ## the xsec used in the templates is FEWZ 3.1
        rescale = self.wxsec(generator)[0]/self.wxsec('fewz3p1')[0]
        values = {}
        histo_file = ROOT.TFile(infile, 'READ')
        for pol in polarizations:
            systs = ['alphaSUp','alphaSDown'] if doAlphaS else [pol+qcdscale+str(ptbin)+charge+direction for qcdscale in ['muR','muF','muRmuF'] for ptbin in xrange(1,11) for direction in ['Up','Down']]
            cp = '{ch}_{pol}'.format(ch=charge,pol=pol)
            xsecs = []
            for iv, val in enumerate(ybins[cp][:-1]):
                if val in excludeYbins: continue
                binname = 'x_W{ch}_{pol}_Ybin_{iy}'.format(ch=charge,pol=pol,iy=iv)
                histo_nominal = histo_file.Get(binname)
                nominal = histo_nominal.Integral()/36000./float(nchannels)
                envelope = 0
                for s in systs:
                    histo_systs = histo_file.Get(binname+'_'+s)
                    envelope = max(envelope, abs(histo_systs.Integral()/36000./float(nchannels) - nominal))
                    if doAlphaS: envelope = 1.5*envelope # one sigma corresponds to +-0.0015 (weights correspond to +-0.001) 
                    envelope = rescale*envelope
                xsecs.append(float(envelope))
            values[pol] = xsecs
        histo_file.Close()
        return values

    def getPDFbandFromXsec(self, histoPDF, charge, infile, netabins, nptbins, firstPtBin=0, histoTotTheory = 0):

        # adding also other alphaS (envelope of Up/Down) and QCD scales (envelope) in quadrature

        print("Inside getPDFbandFromXsec() ...")
        histo_file = ROOT.TFile(infile, 'READ')            

        # NOTE: hnomi is a TH1
        for ieta in range(netabins):
            for ipt in range(firstPtBin,nptbins):

                nomi = "x_W{ch}_lep_ieta_{ie}_ipt_{ip}".format(ch=charge, ie=ieta, ip=ipt)
                hnomi = histo_file.Get(nomi)
                if not hnomi:
                    print("Error in getPDFbandFromXsec(): I couldn't find histogram " + nomi)
                    quit()

                hpdftmp = None            
                pdfquadrsum = 0.0 
                xsecnomi = hnomi.Integral(0, 1+hnomi.GetNbinsX())                

                for ipdf in range(1, 61):
                    pdfvar = nomi + '_pdf{ip}Up'.format(ip=ipdf)
                    hpdftmp = histo_file.Get(pdfvar)
                    if not hpdftmp:
                        print("Error in getPDFbandFromXsec(): I couldn't find histogram " + pdfvar)
                        quit()                    
                    tmpval = hpdftmp.Integral(0, 1+hpdftmp.GetNbinsX()) - xsecnomi                     
                    pdfquadrsum += tmpval * tmpval

                envelopeQCD = 0
                for iqcd in range(1,11):
                    for idir in ["Up", "Down"]:
                        for qcdType in ["muF", "muR", "muRmuF"]:
                            qcdvar = nomi + '_{qt}{iq}{ch}{d}'.format(qt=qcdType,iq=iqcd,ch=charge,d=idir)
                            hqcdtmp = histo_file.Get(qcdvar)
                            if not hqcdtmp:
                                print("Error in getPDFbandFromXsec(): I couldn't find histogram " + qcdvar)
                                quit()                    
                            tmpval = hqcdtmp.Integral(0, 1+hqcdtmp.GetNbinsX()) - xsecnomi                     
                            envelopeQCD = max(abs(envelopeQCD),abs(tmpval))
                                
                alphaS = 0
                for idir in ["Up", "Down"]:
                    alphaSvar = nomi + '_alphaS{d}'.format(d=idir)
                    halphaStmp = histo_file.Get(alphaSvar)
                    if not halphaStmp:
                        print("Error in getPDFbandFromXsec(): I couldn't find histogram " + alphaSdvar)
                        quit()                    
                    # 1.5 is because weight corresponds to 0.001, but should be 0.0015 
                    tmpval = halphaStmp.Integral(0, 1+halphaStmp.GetNbinsX()) - xsecnomi  
                    alphaS = max(abs(alphaS),abs(1.5*tmpval))

                if histoTotTheory:
                    histoPDF.SetBinError(ieta+1, ipt+1, math.sqrt(pdfquadrsum + alphaS * alphaS))
                    histoPDF.SetBinContent(ieta+1, ipt+1, xsecnomi)
                    histoTotTheory.SetBinError(ieta+1, ipt+1, math.sqrt(pdfquadrsum + envelopeQCD * envelopeQCD + alphaS * alphaS))
                    histoTotTheory.SetBinContent(ieta+1, ipt+1, xsecnomi)
                else:
                    histoPDF.SetBinError(ieta+1, ipt+1, math.sqrt(pdfquadrsum + envelopeQCD * envelopeQCD + alphaS * alphaS))
                    histoPDF.SetBinContent(ieta+1, ipt+1, xsecnomi)
                                      
                        
        histo_file.Close()
        return 0


    def getPDFbandFromXsecEta(self, histoPDF, charge, infile, netabins, nptbins, firstPtBin=0):

        # obsolete, use the next one, generic for pt and eta
        print("Inside getPDFbandFromXsecEta() ...")
        histo_file = ROOT.TFile(infile, 'READ')            

        for ieta in range(netabins):
            pdfquadrsum = 0.0
            xsecnomi = 0.0
            xsecpdf = [0.0 for i in range(60)] 
            for ipt in range(firstPtBin,nptbins):

                nomi = "x_W{ch}_lep_ieta_{ie}_ipt_{ip}".format(ch=charge, ie=ieta, ip=ipt)
                hnomi = histo_file.Get(nomi)
                if not hnomi:
                    print("Error in getPDFbandFromXsecEta(): I couldn't find histogram " + nomi)
                    quit()

                hpdftmp = None            
                xsecnomi += hnomi.Integral(0, 1+hnomi.GetNbinsX())                

                for ipdf in range(1, 61):
                    pdfvar = nomi + '_pdf{ip}Up'.format(ip=ipdf)
                    hpdftmp = histo_file.Get(pdfvar)
                    if not hpdftmp:
                        print("Error in getPDFbandFromXsecEta(): I couldn't find histogram " + pdfvar)
                        quit()
                    
                    xsecpdf[ipdf-1] += hpdftmp.Integral(0, 1+hpdftmp.GetNbinsX())

            pdfquadrsum = 0.0
            for ipdf in range(60):
                tmpval = xsecpdf[ipdf] - xsecnomi 
                pdfquadrsum += tmpval*tmpval
            histoPDF.SetBinError(ieta+1, math.sqrt(pdfquadrsum))
            histoPDF.SetBinContent(ieta+1, xsecnomi)                        

        histo_file.Close()
        return 0

    # generalizing previous function getPDFbandFromXsecEta()
    def getPDFbandFromXsec1D(self, histoPDF, charge, infile, netabins, nptbins, firstVarBin=0, firstOtherVarBin=0, isEta = True, histoTotTheory = 0):

        print("Inside getPDFbandFromXsec1D() ...")
        histo_file = ROOT.TFile(infile, 'READ')            

        nvarbins = netabins if isEta else nptbins
        nothervarbins = nptbins if isEta else netabins

        for ivar in range(firstVarBin,nvarbins):
            pdfquadrsum = 0.0
            xsecnomi = 0.0
            xsecpdf = [0.0 for i in range(60)]
            # beofre taking envelope for the following must sum the relevant bins being integrated
            xsecalphaS = {"Up": 0.0, "Down" : 0.0}
            xsecqcd = {}
            for iqcd in range(1,11):
                for idir in ["Up", "Down"]:
                    for itype in ["muR", "muF", "muRmuF"]:
                        xsecqcd["{t}{n}{d}".format(t=itype,n=iqcd,d=idir)] = 0.0
            for iother in range(firstOtherVarBin,nothervarbins):

                nomi = "x_W{ch}_lep_ieta_{ie}_ipt_{ip}".format(ch=charge, ie=ivar if isEta else iother, ip=iother if isEta else ivar)
                hnomi = histo_file.Get(nomi)
                if not hnomi:
                    print("Error in getPDFbandFromXsec1D(): I couldn't find histogram " + nomi)
                    quit()

                hpdftmp = None            
                xsecnomi += hnomi.Integral(0, 1+hnomi.GetNbinsX())                

                for ipdf in range(1, 61):
                    pdfvar = nomi + '_pdf{ip}Up'.format(ip=ipdf)
                    hpdftmp = histo_file.Get(pdfvar)
                    if not hpdftmp:
                        print("Error in getPDFbandFromXsec1D(): I couldn't find histogram " + pdfvar)
                        quit()
                    
                    xsecpdf[ipdf-1] += hpdftmp.Integral(0, 1+hpdftmp.GetNbinsX())

                #########################
                for iqcd in range(1,11):
                    for idir in ["Up", "Down"]:
                        for qcdType in ["muF", "muR", "muRmuF"]:
                            qcdvar = nomi + '_{qt}{iq}{ch}{d}'.format(qt=qcdType,iq=iqcd,ch=charge,d=idir)
                            hqcdtmp = histo_file.Get(qcdvar)
                            if not hqcdtmp:
                                print("Error in getPDFbandFromXsec(): I couldn't find histogram " + qcdvar)
                                quit()                    
                            tmpval = hqcdtmp.Integral(0, 1+hqcdtmp.GetNbinsX())
                            xsecqcd["{qt}{iq}{d}".format(qt=qcdType,iq=iqcd,d=idir)] += tmpval
                                
                for idir in ["Up", "Down"]:
                    alphaSvar = nomi + '_alphaS{d}'.format(d=idir)
                    halphaStmp = histo_file.Get(alphaSvar)
                    if not halphaStmp:
                        print("Error in getPDFbandFromXsec(): I couldn't find histogram " + alphaSdvar)
                        quit()                    
                    tmpval = halphaStmp.Integral(0, 1+halphaStmp.GetNbinsX())
                    xsecalphaS[idir] +=  tmpval
                #########################


            pdfquadrsum = 0.0
            for ipdf in range(60):
                tmpval = xsecpdf[ipdf] - xsecnomi 
                pdfquadrsum += tmpval*tmpval

            envelopeQCD = 0
            for iqcd in range(1,11):
                for idir in ["Up", "Down"]:
                    for qcdType in ["muF", "muR", "muRmuF"]:
                        tmpval = xsecqcd["{qt}{iq}{d}".format(qt=qcdType,iq=iqcd,d=idir)] - xsecnomi
                        envelopeQCD = max(abs(envelopeQCD),abs(tmpval))

            alphaS = 0
            for idir in ["Up", "Down"]:
                # 1.5 is because weight corresponds to 0.001, but should be 0.0015 
                tmpval = xsecalphaS[idir] - xsecnomi
                alphaS = max(abs(alphaS),abs(1.5*tmpval))

            if histoTotTheory:
                histoPDF.SetBinError(ivar+1, math.sqrt(pdfquadrsum + alphaS * alphaS))
                histoPDF.SetBinContent(ivar+1, xsecnomi)
                histoTotTheory.SetBinError(ivar+1, math.sqrt(pdfquadrsum + envelopeQCD * envelopeQCD + alphaS * alphaS))
                histoTotTheory.SetBinContent(ivar+1, xsecnomi)
            else:
                histoPDF.SetBinError(ivar+1, math.sqrt(pdfquadrsum + envelopeQCD * envelopeQCD + alphaS * alphaS))
                histoPDF.SetBinContent(ivar+1, xsecnomi)


        histo_file.Close()
        return 0

#######################

    def getTheoryHistDiffXsecFast(self, xsecWithWptWeights=True, ptmin=-1.0):

        # using native MC@NLO xsec (60400 pb instead of (3*20508.9)pb from fewz3.1)
        # to use fewz3.1, remove "_nativeMCatNLOxsec" from file name

        # adding also other alphaS (envelope of Up/Down) and QCD scales (envelope) in quadrature

        print("Inside getTheoryHistDiffXsecFast() ...")
        # hardcoded for now
        infile = "/afs/cern.ch/work/m/mciprian/public/whelicity_stuff/xsection_genAbsEtaPt_dressed_binningAnalysis_noWpt_yields_nativeMCatNLOxsec.root"
        if xsecWithWptWeights:
            infile = "/afs/cern.ch/work/m/mciprian/public/whelicity_stuff/xsection_genAbsEtaPt_dressed_mu_binningAnalysis_WptWeights_allQCDscales_yields_nativeMCatNLOxsec.root"
        histo_file = ROOT.TFile(infile, 'READ')            

        htheory = {}
        htheory_1Deta = {}
        htheory_1Dpt = {}

        htheory_xsecnorm = {}
        htheory_1Deta_xsecnorm = {}
        htheory_1Dpt_xsecnorm = {}

        htheory_asym = {}
        htheory_1Deta_asym = {}
        htheory_1Dpt_asym = {}

        print("ABSOLUTE XSEC")
        for charge in ["plus", "minus"]:
 
            nomi = "gen_ptl1_absetal1_dressed_binAna_W{ch}_mu".format(ch=charge)
            hnomi = histo_file.Get(nomi + "_central")
            self.checkHistInFile(hnomi, nomi + "_central", infile, message="in getTheoryHistDiffXsecFast()")
            hnomi.SetDirectory(0)

            for ipdf in range(1,61):
                name = "pdf{ip}".format(ip=ipdf)
                pdfvar = nomi + '_{n}'.format(n=name)
                htheory[(charge,name)] = histo_file.Get(pdfvar)
                self.checkHistInFile(htheory[(charge,name)], pdfvar, infile, message="in getTheoryHistDiffXsecFast()")
                htheory[(charge,name)].SetDirectory(0)

            for idir in ["Up", "Dn"]:  # should have been Down, but the histogram migt have Dn
                for qcdType in ["muF", "muR", "muRmuF"]:
                    name = "{qt}{d}".format(qt=qcdType,d=idir)
                    qcdvar = nomi + '_{n}'.format(n=name)
                    htheory[(charge,name)] = histo_file.Get(qcdvar)
                    self.checkHistInFile(htheory[(charge,name)], qcdvar, infile, message="in getTheoryHistDiffXsecFast()")
                    htheory[(charge,name)].SetDirectory(0)

                name = "alphaS{d}".format(d=idir)
                alphaSvar = nomi + '_{n}'.format(n=name)
                htheory[(charge,name)] = histo_file.Get(alphaSvar)
                self.checkHistInFile(htheory[(charge,name)], alphaSvar, infile, message="in getTheoryHistDiffXsecFast()")
                htheory[(charge,name)].SetDirectory(0)

            htheory[(charge,"nominal")] = hnomi

        # key has both plus and minus charges
        print("NORMALIZED XSEC")
        for key in htheory:
            # this uses only the xsec in acceptance, |eta| in [0, 2.4], pT in [26, 56]
            # WARNING: for electrons I should restrict pt range to [30, 56], unless the fit uses also bkg processes to make the total xsec
            ptminbin = 1
            if float(ptmin) > 0:
                ptminbin = htheory[key].GetYaxis().FindFixBin(float(ptmin) + 0.0001) # adding epsilon for safety
            # projections
            htheory_1Deta[key] = htheory[key].ProjectionX(htheory[key].GetName()+"_1Deta",ptminbin,htheory[key].GetNbinsY())
            htheory_1Deta[key].SetDirectory(0)
            htheory_1Dpt[key] = htheory[key].ProjectionY(htheory[key].GetName()+"_1Dpt",1,htheory[key].GetNbinsX())             
            htheory_1Dpt[key].SetDirectory(0)
            htheory_xsecnorm[key] = htheory[key].Clone(htheory[key].GetName() + "_xsecnorm")
            htheory_xsecnorm[key].SetDirectory(0)
            htheory_xsecnorm[key].Scale(1./htheory_xsecnorm[key].Integral(1, htheory_xsecnorm[key].GetNbinsX(), 
                                                                          ptminbin, htheory_xsecnorm[key].GetNbinsY())) 
            htheory_1Deta_xsecnorm[key] = htheory_1Deta[key].Clone(htheory_1Deta[key].GetName() + "_xsecnorm")
            htheory_1Deta_xsecnorm[key].SetDirectory(0)
            htheory_1Deta_xsecnorm[key].Scale(1./htheory_1Deta_xsecnorm[key].Integral()) 
            htheory_1Dpt_xsecnorm[key] = htheory_1Dpt[key].Clone(htheory_1Dpt[key].GetName() + "_xsecnorm")
            htheory_1Dpt_xsecnorm[key].SetDirectory(0)
            htheory_1Dpt_xsecnorm[key].Scale(1./htheory_1Dpt_xsecnorm[key].Integral(ptminbin,htheory_1Dpt_xsecnorm[key].GetNbinsX())) 

        # now charge asymmetry
        print("CHARGE AYMMETRY")
        htmp_num = htheory[("plus","nominal")].Clone("_tmp_helper_asym_num")
        htmp_den = htheory[("plus","nominal")].Clone("_tmp_helper_asym_den")
        htmp_num_1Deta = htheory_1Deta[("plus","nominal")].Clone("_tmp_helper_asym_num_1Deta")
        htmp_den_1Deta = htheory_1Deta[("plus","nominal")].Clone("_tmp_helper_asym_den_1Deta")
        htmp_num_1Dpt = htheory_1Dpt[("plus","nominal")].Clone("_tmp_helper_asym_num_1Dpt")
        htmp_den_1Dpt = htheory_1Dpt[("plus","nominal")].Clone("_tmp_helper_asym_den_1Dpt")
        for key in htheory:
            # filter keys for plus charge
            if str(key[0]) == "minus": continue
            # 2D
            keyAsym = ("all",key[1]) # to allow for homogeneous treatment as for the charged xsecs
            htheory_asym[keyAsym] = htheory[key].Clone(htheory[key].GetName() + "_asymm")
            htheory_asym[keyAsym].SetDirectory(0)
            htmp_num.Add(htheory[key], htheory[("minus",key[1])], 1.0, -1.0)
            htmp_den.Add(htheory[key], htheory[("minus",key[1])], 1.0,  1.0)
            htheory_asym[keyAsym].Divide(htmp_num,htmp_den)
            # 1D eta
            htheory_1Deta_asym[keyAsym] = htheory_1Deta[key].Clone(htheory_1Deta[key].GetName() + "_asymm")
            htheory_1Deta_asym[keyAsym].SetDirectory(0)
            htmp_num_1Deta.Add(htheory_1Deta[key], htheory_1Deta[("minus",key[1])], 1.0, -1.0)
            htmp_den_1Deta.Add(htheory_1Deta[key], htheory_1Deta[("minus",key[1])], 1.0,  1.0)
            htheory_1Deta_asym[keyAsym].Divide(htmp_num_1Deta,htmp_den_1Deta)
            # 1D pt
            htheory_1Dpt_asym[keyAsym] = htheory_1Dpt[key].Clone(htheory_1Dpt[key].GetName() + "_asymm")
            htheory_1Dpt_asym[keyAsym].SetDirectory(0)
            htmp_num_1Dpt.Add(htheory_1Dpt[key], htheory_1Dpt[("minus",key[1])], 1.0, -1.0)
            htmp_den_1Dpt.Add(htheory_1Dpt[key], htheory_1Dpt[("minus",key[1])], 1.0,  1.0)
            htheory_1Dpt_asym[keyAsym].Divide(htmp_num_1Dpt,htmp_den_1Dpt)

        histo_file.Close()
        ret = { "xsec"     : [htheory, htheory_1Deta, htheory_1Dpt],
                "xsecnorm" : [htheory_xsecnorm, htheory_1Deta_xsecnorm, htheory_1Dpt_xsecnorm],
                "asym"     : [htheory_asym, htheory_1Deta_asym, htheory_1Dpt_asym],
                "listkeys" : [key[1] for key in htheory]}
        return ret

#######################

    def getTheoryBandDiffXsec(self, hretTotTheory, hretPDF, theovars, hists, netabins, nptbins, charge="all"):
        
        # AlphaS and QCD can have asymmetric uncertainties, so I need a TGraphAsymmErrors()

        # charge == all for asymmetry, plus or minus otherwise
        hnomi = hists[(charge, "nominal")]

        for ieta in range(netabins):
            for ipt in range(nptbins):
                pdfquadrsum = 0.0
                envelopeQCD = 0.0     # muR, muF, muRmuF, with Up and Down
                envelopeAlphaS = 0.0  # basically Up and Down
                nomi = hnomi.GetBinContent(ieta+1, ipt+1)
                #envelopeAlphaS = 0.0
                #envelopeQCD = 0.0                
                envelopeQCD_vals = [nomi]     # muR, muF, muRmuF, with Up and Down
                envelopeAlphaS_vals = [nomi]  # basically Up and Down
                for nuis in theovars:
                    if nuis == "nominal": continue
                    if "pdf" in nuis:
                        tmpvar = hists[(charge,nuis)].GetBinContent(ieta+1,ipt+1) - nomi
                        pdfquadrsum += tmpvar * tmpvar
                    elif "alpha" in nuis:
                        # 1.5 is because weight corresponds to 0.001, but should be 0.0015
                        #tmpvar = hists[(charge,nuis)].GetBinContent(ieta+1,ipt+1) - nomi
                        #envelopeAlphaS = max(abs(envelopeAlphaS), 1.5* abs(tmpvar))
                        alt = hists[(charge,nuis)].GetBinContent(ieta+1, ipt+1)
                        altscaled = nomi * math.exp( 1.5 * math.log(alt/nomi) )
                        envelopeAlphaS_vals.append(altscaled)
                    else:
                        #tmpvar = hists[(charge,nuis)].GetBinContent(ieta+1,ipt+1) - nomi
                        #envelopeQCD = max(abs(envelopeQCD), abs(tmpvar))
                        envelopeQCD_vals.append(hists[(charge,nuis)].GetBinContent(ieta+1, ipt+1))

                # pdfAndAlphaQuadSum = pdfquadrsum + envelopeAlphaS * envelopeAlphaS
                # hretPDF.SetBinContent(ieta+1,ipt+1, hnomi.GetBinContent(ieta+1,ipt+1))
                # hretTotTheory.SetBinContent(ieta+1,ipt+1, hnomi.GetBinContent(ieta+1,ipt+1))
                # hretPDF.SetBinError(ieta+1,ipt+1, math.sqrt(pdfAndAlphaQuadSum))
                # hretTotTheory.SetBinError(ieta+1,ipt+1, math.sqrt(pdfAndAlphaQuadSum + envelopeQCD * envelopeQCD))

                # this graph will be used for the unrolled xsec, which is a TH1 drawn with bin width equal to 1
                errX = 0.5
                alphaSErrorHigh = abs(max(envelopeAlphaS_vals)-nomi)
                alphaSErrorLow  = abs(min(envelopeAlphaS_vals)-nomi)
                QCDErrorHigh    = abs(max(envelopeQCD_vals)-nomi)
                QCDErrorLow     = abs(min(envelopeQCD_vals)-nomi)
                pdfAlphaSErrorHigh = math.sqrt(pdfquadrsum + alphaSErrorHigh * alphaSErrorHigh) 
                pdfAlphaSErrorLow  = math.sqrt(pdfquadrsum + alphaSErrorLow * alphaSErrorLow) 
                totTheoryErrorHigh = math.sqrt(pdfquadrsum + alphaSErrorHigh * alphaSErrorHigh + QCDErrorHigh * QCDErrorHigh)
                totTheoryErrorLow  = math.sqrt(pdfquadrsum + alphaSErrorLow * alphaSErrorLow + QCDErrorLow * QCDErrorLow)

                ibin = ieta + ipt * netabins # from 0 to neta*npt-1
                hretPDF.SetPoint(ibin, ibin+1, nomi)
                hretPDF.SetPointError(ibin, errX, errX, pdfAlphaSErrorLow,  pdfAlphaSErrorHigh)
                hretTotTheory.SetPoint(ibin, ibin+1, nomi)
                hretTotTheory.SetPointError(ibin, errX, errX, totTheoryErrorLow, totTheoryErrorHigh)



    def getTheoryBandDiffXsec1Dproj(self, hretTotTheory, hretPDF, theovars, hists, nvarbins, charge="all" ):
        
        # AlphaS and QCD can have asymmetric uncertainties, so I need a TGraphAsymmErrors()

        hnomi = hists[(charge,"nominal")]
        for ivar in range(nvarbins):
            pdfquadrsum = 0.0
            nomi = hnomi.GetBinContent(ivar+1)
            envelopeQCD_vals = [nomi]     # muR, muF, muRmuF, with Up and Down
            envelopeAlphaS_vals = [nomi]  # basically Up and Down
            for nuis in theovars:
                if nuis == "nominal": continue
                if "pdf" in nuis:
                    tmpvar = hists[(charge,nuis)].GetBinContent(ivar+1) - nomi
                    pdfquadrsum += tmpvar * tmpvar
                elif "alpha" in nuis:
                    # 1.5 is because weight corresponds to 0.001, but should be 0.0015
                    alt = hists[(charge,nuis)].GetBinContent(ivar+1)
                    altscaled = nomi * math.exp( 1.5 * math.log(alt/nomi) )
                    envelopeAlphaS_vals.append(altscaled)
                else:
                    envelopeQCD_vals.append(hists[(charge,nuis)].GetBinContent(ivar+1))

            errX = 0.5 * hnomi.GetBinWidth(ivar+1) # symmetric bins in x, error is just half the bin width
            alphaSErrorHigh = abs(max(envelopeAlphaS_vals)-nomi)
            alphaSErrorLow  = abs(min(envelopeAlphaS_vals)-nomi)
            QCDErrorHigh    = abs(max(envelopeQCD_vals)-nomi)
            QCDErrorLow     = abs(min(envelopeQCD_vals)-nomi)
            pdfAlphaSErrorHigh = math.sqrt(pdfquadrsum + alphaSErrorHigh * alphaSErrorHigh) 
            pdfAlphaSErrorLow  = math.sqrt(pdfquadrsum + alphaSErrorLow * alphaSErrorLow) 
            totTheoryErrorHigh = math.sqrt(pdfquadrsum + alphaSErrorHigh * alphaSErrorHigh + QCDErrorHigh * QCDErrorHigh)
            totTheoryErrorLow  = math.sqrt(pdfquadrsum + alphaSErrorLow * alphaSErrorLow + QCDErrorLow * QCDErrorLow)

            hretPDF.SetPoint(ivar, hnomi.GetBinCenter(ivar+1), nomi)
            hretPDF.SetPointError(ivar, errX, errX, pdfAlphaSErrorLow,  pdfAlphaSErrorHigh)
            hretTotTheory.SetPoint(ivar, hnomi.GetBinCenter(ivar+1), nomi)
            hretTotTheory.SetPointError(ivar, errX, errX, totTheoryErrorLow, totTheoryErrorHigh)


    def checkTheoryBandDiffXsec1Dproj(self, hretPDF, hretAlpha, hretQCD, theovars, hists, nvarbins, charge="all" ):
        
        # hretAlpha, hretQCD can have asymmetric uncertainties, so I need a TGraphAsymmErrors()
        # pdfs can stay as a TH1

        # to test width of the bands
        hnomi = hists[(charge,"nominal")]
        for ivar in range(nvarbins):
            pdfquadrsum = 0.0
            nomi = hnomi.GetBinContent(ivar+1)
            envelopeQCD_vals = [nomi]     # muR, muF, muRmuF, with Up and Down
            envelopeAlphaS_vals = [nomi]  # basically Up and Down
            for nuis in theovars:
                if nuis == "nominal": continue
                if "pdf" in nuis:
                    tmpvar = hists[(charge,nuis)].GetBinContent(ivar+1) - nomi
                    pdfquadrsum += tmpvar * tmpvar
                elif "alpha" in nuis:
                    # 1.5 is because weight corresponds to 0.001, but should be 0.0015
                    #tmpvar = hists[(charge,nuis)].GetBinContent(ivar+1) - hnomi.GetBinContent(ivar+1)
                    #envelopeAlphaS = max(abs(envelopeAlphaS), 1.5* abs(tmpvar))
                    # if alphaHist is larger than nomi difference is positive and I add the variations time 1.5, otherwise
                    # it is negative and I subtract
                    alt = hists[(charge,nuis)].GetBinContent(ivar+1)
                    altscaled = nomi * math.exp( 1.5 * math.log(alt/nomi) )
                    envelopeAlphaS_vals.append(altscaled)
                else:
                    #tmpvar = hists[(charge,nuis)].GetBinContent(ivar+1) - hnomi.GetBinContent(ivar+1)
                    #envelopeQCD = max(abs(envelopeQCD), abs(tmpvar))
                    envelopeQCD_vals.append(hists[(charge,nuis)].GetBinContent(ivar+1))

            hretPDF.SetBinContent(ivar+1, nomi)
            hretPDF.SetBinError(ivar+1, math.sqrt(pdfquadrsum))
            errX = 0.5 * hnomi.GetBinWidth(ivar+1) # symmetric bins in x, error is just half the bin width
            hretAlpha.SetPoint(ivar, hnomi.GetBinCenter(ivar+1), nomi)
            hretAlpha.SetPointError(ivar, errX, errX, abs(min(envelopeAlphaS_vals)-nomi), abs(max(envelopeAlphaS_vals)-nomi))
            hretQCD.SetPoint(ivar, hnomi.GetBinCenter(ivar+1), nomi)
            hretQCD.SetPointError(ivar, errX, errX, abs(min(envelopeQCD_vals)-nomi), abs(max(envelopeQCD_vals)-nomi))
            #hretAlpha.SetBinContent(ivar+1, hnomi.GetBinContent(ivar+1))
            #hretAlpha.SetBinError(ivar+1, envelopeAlphaS)
            #hretQCD.SetBinContent(ivar+1, hnomi.GetBinContent(ivar+1))
            #hretQCD.SetBinError(ivar+1, envelopeQCD)


#######################

    def getParametersFromWS(self, ws, regexp):

        ## get all the nuisance parameters from the workspace
        pars = ws.allVars()
        pars = ROOT.RooArgList(pars)
        ## this has to be a loop over a range... doesn't work otherwise
        parameters = []
        all_parameters = []
        pois_regexps = list(regexp.split(','))
        ## get the parameters to scan from the list of allVars and match them
        ## to the given regexp
        for i in range(len(pars)):
            tmp_name = pars[i].GetName()
            if '_In' in tmp_name: ## those are the input parameters
                continue
            if tmp_name in ['CMS_th1x', 'r']: ## don't want those
                continue
            all_parameters.append(tmp_name)
            for poi in pois_regexps:
                if re.match(poi, tmp_name):
                    parameters.append(pars[i].GetName())

        return parameters

    def translateWStoTF(self, pname):
        if not 'r_' in pname:
            return pname
     
        if 'Wplus' in pname or 'Wminus' in pname:
            pnew = pname.replace('r_','')
            pnew = pnew.split('_')
            pnew.insert(-2, 'mu' )#if 'Wmu' in options.tensorflow else 'el')
            pnew = '_'.join(pnew)
            return pnew
     
        return -1

    def getFromHessian(self, infile, keepGen=False, takeEntry=0, params=[]):
        _dict = {}
        
        f = ROOT.TFile(infile, 'read')
        tree = f.Get('fitresults')
        if tree.GetEntries() > 1:
            print('YOUR INPUT FILE HAS MORE THAN ONE FIT INSIDE. THIS IS PROBABLY NOT A HESSIAN FILE!!!')
            print('will take the first entry by default. unless specified otherwise')
            ##sys.exit()
        lok  = tree.GetListOfLeaves()
        for p in lok:
            if '_err'   in p.GetName(): continue
            if '_minos' in p.GetName(): continue
            if '_gen'   in p.GetName() and not keepGen: continue
            if '_In'    in p.GetName(): continue

            if len(params):
                match = [re.match(param,p.GetName()) for param in params]
                if not any(match): continue

            if not p.GetName()+'_err' in lok and not keepGen: continue

            for iev,ev in enumerate(tree):
                if not iev == takeEntry: continue ## this should work... i guess.
                mean = getattr(ev, p.GetName())
                err  = getattr(ev, p.GetName()+'_err') if hasattr(ev, p.GetName()+'_err') else 0

            _dict[p.GetName()] = (mean, mean+err, mean-err)
     
        return _dict


    def getFromToys(self, infile, keepGen=False, params=[]):
        _dict = {}
        
        f = ROOT.TFile(infile, 'read')
        tree = f.Get('fitresults')
        lok  = tree.GetListOfLeaves()
        
        for p in lok:
            if '_err'   in p.GetName(): continue
            if '_minos' in p.GetName(): continue
            if '_gen'   in p.GetName() and not keepGen: continue
            if '_In'    in p.GetName(): continue

            if len(params):
                match = [re.match(param,p.GetName()) for param in params]
                if not any(match): continue

            print('gettin parameter ', p.GetName(), 'from toys file')
            
            tmp_hist = ROOT.TH1F(p.GetName(),p.GetName(), 100000, -5000., 5000.)
            tree.Draw(p.GetName()+'>>'+p.GetName())
            #tmp_hist = ROOT.gPad.GetPrimitive('foob')
            mean = tmp_hist.GetMean()
            err  = tmp_hist.GetRMS()
            _dict[p.GetName()] = (mean, mean+err, mean-err)
            del tmp_hist

        return _dict

    def getHistosFromToys(self, infile, nbins=100, xlow=-3.0, xup=3.0, getPull=False, matchBranch=None,excludeBranch=None, selection="", 
                          setStatOverflow=False, getMedian=False):

        # getPull = True will return a histogram centered at 0 and with expected rms=1, obtained as (x-x_gen)/x_err

        _dict = {}
        
        f = ROOT.TFile(infile, 'read')
        tree = f.Get('fitresults')
        lok  = tree.GetListOfLeaves()

        #nMaxBranch = 10
        #np = 0
        im = 1  # for median

        tree.SetBranchStatus("*",0)  # disabling and enabling branches makes the loop on events faster, at least if the median is used

        for p in lok:

            tree.SetBranchStatus(p.GetName(),1)
            
            if '_err'   in p.GetName(): continue
            if '_minos' in p.GetName(): continue
            if '_gen'   in p.GetName(): continue
            if '_In'    in p.GetName(): continue
            if matchBranch and not any(re.match(poi,p.GetName()) for poi in matchBranch.split(',')): continue
            if excludeBranch and any(re.match(excl,p.GetName()) for excl in excludeBranch.split(',')): continue

            # mainly for tests
            #if np == nMaxBranch: break
            #np += 1
            
            #print("Loading parameter --> %s " % p.GetName())
            #self.cmssw_version = os.environ['CMSSW_VERSION']
            #self.isRecentRelease = (len(self.cmssw_version) and int(self.cmssw_version.split('_')[1]) > 8)

            if getPull and (p.GetName()+"_gen") in lok and (p.GetName()+"_err") in lok:                
                tree.SetBranchStatus(p.GetName()+"_gen",1)
                tree.SetBranchStatus(p.GetName()+"_err",1)
                #print(" Making pull --> (x-x_gen)/x_err for parameter %s" % p.GetName())
                tmp_hist_tmp = ROOT.TH1F(p.GetName()+"_tmp",p.GetName()+"_tmp", nbins, xlow, xup)
                if self.isRecentRelease:
                    if setStatOverflow: tmp_hist_tmp.SetStatOverflows(1)
                    else              : tmp_hist_tmp.SetStatOverflows(0)
                tmp_hist = ROOT.TH1F(p.GetName(),p.GetName(), 100, -3, 3)
                expression = "({p}-{pgen})/{perr}".format(p=p.GetName(),pgen=p.GetName()+"_gen",perr=p.GetName()+"_err")
                tree.Draw(expression+'>>'+p.GetName(),selection)
                tree.Draw(p.GetName()+'>>'+p.GetName()+"_tmp",selection)
                mean = tmp_hist_tmp.GetMean()
                err  = tmp_hist_tmp.GetRMS()
            else:
                tmp_hist = ROOT.TH1F(p.GetName(),p.GetName(), nbins, xlow, xup)
                if self.isRecentRelease:
                    if setStatOverflow: tmp_hist.SetStatOverflows(1)
                    else              : tmp_hist.SetStatOverflows(0)
                tree.Draw(p.GetName()+'>>'+p.GetName(),selection)
                mean = tmp_hist.GetMean()
                err  = tmp_hist.GetRMS()
            tmp_hist.SetDirectory(None)

            if getMedian:
                print("{n}) Computing median for {pn}".format(n=im,pn=p.GetName()))
                im += 1
                vals = []
                #binCount = 1
                tot = tree.GetEntries()
                for ev in tree:
                    #sys.stdout.write('Bin {num}/{tot}   \r'.format(num=binCount,tot=tot))
                    #sys.stdout.flush()
                    #binCount += 1
                    vals.append(getattr(ev, p.GetName()))
                    
                vals.sort()
                nElements = len(vals)            
                if nElements%2: median = vals[(nElements-1)/2]
                else:           median = 0.5 * (vals[nElements/2] +  vals[nElements/2 + 1])
                _dict[p.GetName()] = (median, median+err, median-err, tmp_hist)
            else:
                _dict[p.GetName()] = (mean, mean+err, mean-err, tmp_hist)

            tree.SetBranchStatus(p.GetName(),0)
            tree.SetBranchStatus(p.GetName()+"_gen",0)
            tree.SetBranchStatus(p.GetName()+"_err",0)
     
        return _dict



    def getExprFromToys(self, name, expression, infile):
        f = ROOT.TFile(infile, 'read')
        tree = f.Get('fitresults')        
        tmp_hist = ROOT.TH1F(name,name, 100000, -100., 5000.)
        tree.Draw(expression+'>>'+name)
        mean = tmp_hist.GetMean()
        err  = tmp_hist.GetRMS()
        return (mean, mean+err, mean-err)


    def getNormalizedXsecFromToys(self, ybins, charge, pol, channel, iy, infile, absYmax=6.0):
        cp = '{ch}_{pol}'.format(ch=charge,pol=pol)
        ybins_expr = []
        for allpol in ['left','right', 'long']:
            for iv, val in enumerate(ybins[cp][:-1]):
                if abs(val)<absYmax:
                    ybins_expr.append('W{charge}_{pol}_W{charge}_{pol}_{ch}_Ybin_{iy}_pmaskedexp'.format(charge=charge,pol=allpol,ch=channel,iy=iv))
        num = 'W{charge}_{pol}_W{charge}_{pol}_{ch}_Ybin_{iy}_pmaskedexp'.format(charge=charge,pol=pol,ch=channel,iy=iy)
        den = '('+'+'.join(ybins_expr)+')'        
        ret = self.getExprFromToys(charge+pol+channel+str(iy),'{num}/{den}'.format(num=num,den=den),infile)
        return ret

    def getAsymmetryFromToys(self, pol, channel, iy, infile):
        expr = '(Wplus_{pol}_Wplus_{pol}_{ch}_Ybin_{iy}_pmaskedexp - Wminus_{pol}_Wminus_{pol}_{ch}_Ybin_{iy}_pmaskedexp)/(Wplus_{pol}_Wplus_{pol}_{ch}_Ybin_{iy}_pmaskedexp + Wminus_{pol}_Wminus_{pol}_{ch}_Ybin_{iy}_pmaskedexp)'.format(pol=pol,ch=channel,iy=iy)
        ret = self.getExprFromToys('chargeAsym',expr,infile)
        return ret

    def getDiffXsecAsymmetryFromToys(self, ieta, ipt, infile):
        xplus  = "Wplus_lep_ieta_{ieta}_ipt_{ipt}_pmaskedexp".format(ieta=ieta,ipt=ipt)
        xminus = "Wminus_lep_ieta_{ieta}_ipt_{ipt}_pmaskedexp".format(ieta=ieta,ipt=ipt)
        expr = '({pl}-{mn})/({pl}+{mn})'.format(pl=xplus,mn=xminus)
        ret = self.getExprFromToys('chargeAsym',expr,infile)
        return ret

    def getDenExpressionForNormDiffXsec(self, charge, netabins, nptbins):
        binsToNormalize = []
        for ieta in range(netabins):
            for ipt in range(nptbins):
                denChunk = "W{c}_lep_ieta_{ieta}_ipt_{ipt}_pmaskedexp".format(c=charge,ieta=ieta,ipt=ipt)
                binsToNormalize.append(denChunk)
        den = "+".join(x for x in binsToNormalize)
        den = "(" + den + ")"
        return den

    def getNormalizedDiffXsecFromToys(self, charge, ieta, ipt, infile, den, friendTree=""):
        num = "W{c}_lep_ieta_{ieta}_ipt_{ipt}_pmaskedexp".format(c=charge,ieta=ieta,ipt=ipt)
        expr = '{num}/{den}'.format(num=num,den=den)
        ret = self.getExprFromToys('normDiffXsec',expr,infile)
        return ret


    def getDiffXsecFromToys(self, charge, ieta, ipt, infile):
        expr = "W{c}_lep_ieta_{ieta}_ipt_{ipt}_pmaskedexp".format(c=charge,ieta=ieta,ipt=ipt)
        ret = self.getExprFromToys('diffXsec',expr,infile)
        return ret

    def getExprFromToysFast(self, name, expression, nHistBins=100000, minHist=-100., maxHist=5000., tree=None):
        tmp_hist = ROOT.TH1F(name,name, nHistBins, minHist, maxHist)
        #self.cmssw_version = os.environ['CMSSW_VERSION']
        #self.isRecentRelease = (len(self.cmssw_version) and int(self.cmssw_version.split('_')[1]) > 8)
        if self.isRecentRelease: tmp_hist.SetStatOverflows(1)
        tree.Draw(expression+'>>'+name)
        mean = tmp_hist.GetMean()
        err  = tmp_hist.GetRMS()
        return (mean, mean+err, mean-err)

    def getNormalizedDiffXsecFromToysFast(self, charge, ieta, ipt, den, nHistBins=1000, minHist=0., maxHist=0.1, tree=None):
        num = "W{c}_lep_ieta_{ieta}_ipt_{ipt}_pmaskedexp".format(c=charge,ieta=ieta,ipt=ipt)
        expr = '{num}/{den}'.format(num=num,den=den)
        ret = self.getExprFromToysFast('normDiffXsec',expr, nHistBins, minHist, maxHist, tree=tree)
        return ret

    def getDiffXsecFromToysFast(self, charge, ieta, ipt, nHistBins=2000, minHist=0., maxHist=200., tree=None):
        expr = "W{c}_lep_ieta_{ieta}_ipt_{ipt}_pmaskedexp".format(c=charge,ieta=ieta,ipt=ipt)
        ret = self.getExprFromToysFast('diffXsec',expr,nHistBins, minHist, maxHist, tree=tree)
        return ret

    def getSignalStrengthFromToysFast(self, channel, charge, ieta, ipt, nHistBins=400, minHist=0., maxHist=2., tree=None):
        expr = "W{c}_lep_ieta_{ieta}_ipt_{ipt}_mu".format(c=charge,ieta=ieta,ipt=ipt)
        ret = self.getExprFromToysFast('diffXsec',expr, nHistBins, minHist, maxHist, tree=tree)
        return ret


    def getDiffXsecAsymmetryFromToysFast(self, ieta, ipt, nHistBins=2000, minHist=0., maxHist=1.0, tree=None):
        xplus  = "Wplus_lep_ieta_{ieta}_ipt_{ipt}_pmaskedexp".format(ieta=ieta,ipt=ipt)
        xminus = "Wminus_lep_ieta_{ieta}_ipt_{ipt}_pmaskedexp".format(ieta=ieta,ipt=ipt)
        expr = '({pl}-{mn})/({pl}+{mn})'.format(pl=xplus,mn=xminus)
        ret = self.getExprFromToysFast('chargeAsym',expr, nHistBins, minHist, maxHist, tree=tree)
        return ret

    def getExprFromHessian(self, name, expression, infile):
        f = ROOT.TFile(infile, 'read')
        tree = f.Get('fitresults')        
        tmp_hist = ROOT.TH1F(name,name, 100000, -100., 5000.)
        tree.Draw(expression+'>>'+name)
        mean = tmp_hist.GetMean()  # if this is hessian and not toys, there is just one entry, so the mean is the entry
        #err  = tmp_hist.GetRMS()  # not used in this context, (we are going to use this expression mainly for charge asymmetry, the uncertainty must be taken from toys)
        return mean

    def getDiffXsecAsymmetryFromHessian(self, ieta, ipt, infile):
        xplus  = "Wplus_lep_ieta_{ieta}_ipt_{ipt}_pmaskedexp".format(ieta=ieta,ipt=ipt)
        xminus = "Wminus_lep_ieta_{ieta}_ipt_{ipt}_pmaskedexp".format(ieta=ieta,ipt=ipt)
        expr = '({pl}-{mn})/({pl}+{mn})'.format(pl=xplus,mn=xminus)
        ret = self.getExprFromHessian('chargeAsym',expr,infile)
        return ret


    def getNormalizedDiffXsecFromHessian(self, charge, ieta, ipt, infile,den):
        num = "W{c}_lep_ieta_{ieta}_ipt_{ipt}_pmaskedexp".format(c=charge,ieta=ieta,ipt=ipt)
        expr = '{num}/{den}'.format(num=num,den=den)
        ret = self.getExprFromHessian('normDiffXsec',expr,infile)
        return ret

    def getDiffXsecFromHessian(self, charge, ieta, ipt, infile):
        expr = "W{c}_lep_ieta_{ieta}_ipt_{ipt}_pmaskedexp".format(c=charge,ieta=ieta,ipt=ipt)
        ret = self.getExprFromHessian('diffXsec',expr,infile)
        return ret

    ######## HESSIAN FAST

    def getExprFromHessianFast(self, name, expression, nHistBins=100000, minHist=-100., maxHist=5000., tree=None):
        tmp_hist = ROOT.TH1F(name,name, nHistBins, minHist, maxHist)
        if self.isRecentRelease: tmp_hist.SetStatOverflows(1)
        tree.Draw(expression+'>>'+name)
        mean = tmp_hist.GetMean()  # if this is hessian and not toys, there is just one entry, so the mean is the entry
        return mean

    def getDiffXsecAsymmetryFromHessianFast(self, ieta, ipt, nHistBins=2000, minHist=0., maxHist=1.0, tree=None, getErr=False, getGen=False):
        expr = "W_lep_ieta_{ieta}_ipt_{ipt}_chargeasym".format(ieta=ieta,ipt=ipt)
        if getErr: expr += "_err"
        elif getGen: expr += "_gen"
        ret = self.getExprFromHessianFast('chargeAsym',expr, nHistBins, minHist, maxHist, tree=tree)
        return ret


    def getNormalizedDiffXsecFromHessianFast(self, charge, ieta, ipt, nHistBins=1000, minHist=0., maxHist=0.1, tree=None, getErr=False, getGen=False):
        expr = "W{c}_lep_ieta_{ieta}_ipt_{ipt}_pmaskedexpnorm".format(c=charge,ieta=ieta,ipt=ipt)
        if getErr: expr += "_err"
        elif getGen: expr += "_gen"
        ret = self.getExprFromHessianFast('normDiffXsec',expr, nHistBins, minHist, maxHist, tree=tree)
        return ret

    def getDiffXsecFromHessianFast(self, charge, ieta, ipt,  nHistBins=2000, minHist=0., maxHist=200., tree=None, getErr=False, getGen=False):
        expr = "W{c}_lep_ieta_{ieta}_ipt_{ipt}_pmaskedexp".format(c=charge,ieta=ieta,ipt=ipt)
        if getErr: expr += "_err"
        elif getGen: expr += "_gen"
        ret = self.getExprFromHessianFast('diffXsec',expr, nHistBins, minHist, maxHist, tree=tree)
        return ret


    # this is for single charge
    def getDiffXsec1DFromHessianFast(self, charge, ietaORipt, isIeta=True, nHistBins=5000, minHist=0., maxHist=5000., tree=None, getErr=False, getGen=False):
        expr = "W{c}_lep_i{var}_{ivar}_sumxsec".format(c=charge,var="eta" if isIeta else "pt", ivar=ietaORipt)
        if getErr: expr += "_err"
        elif getGen: expr += "_gen"
        ret = self.getExprFromHessianFast('diffXsec1D_{var}'.format(var="eta" if isIeta else "pt"),expr, nHistBins, minHist, maxHist, tree=tree)
        return ret

    def getNormalizedDiffXsec1DFromHessianFast(self, charge, ietaORipt, isIeta=True, nHistBins=1000, minHist=0., maxHist=1., tree=None, getErr=False, getGen=False):
        expr = "W{c}_lep_i{var}_{ivar}_sumxsecnorm".format(c=charge,var="eta" if isIeta else "pt", ivar=ietaORipt)
        if getErr: expr += "_err"
        elif getGen: expr += "_gen"
        ret = self.getExprFromHessianFast('diffXsecNorm1D_{var}'.format(var="eta" if isIeta else "pt"),expr, nHistBins, minHist, maxHist, tree=tree)
        return ret


    def getDiffXsecAsymmetry1DFromHessianFast(self, ietaORipt, isIeta=True, 
                                              nHistBins=1000, minHist=0., maxHist=1., tree=None, getErr=False, getGen=False):
        expr = "W_lep_i{var}_{ivar}_chargemetaasym".format(var="eta" if isIeta else "pt", ivar=ietaORipt)
        if getErr: expr += "_err"
        elif getGen: expr += "_gen"
        ret = self.getExprFromHessianFast('chargeAsym1D_{var}'.format(var="eta" if isIeta else "pt"),expr, nHistBins, minHist, maxHist, tree=tree)
        return ret


    def getSignalStrengthFromHessianFast(self, charge, ieta, ipt, nHistBins=200, minHist=0.5, maxHist=1.5, tree=None, getErr=False, getGen=False):
        expr = "W{c}_lep_ieta_{ieta}_ipt_{ipt}_mu".format(c=charge,ieta=ieta,ipt=ipt)
        if getErr: expr += "_err"
        elif getGen: expr += "_gen"
        ret = self.getExprFromHessianFast('normDiffXsec',expr, nHistBins, minHist, maxHist, tree=tree)
        return ret


    #################################### 

    def getFromScans(self, indir):
        _dict = {}
        
        for sd in os.listdir(indir):
     
            if 'jobs' in sd: continue
     
            par = self.translateWStoTF(sd) ## parameter name different than in TF
            f = ROOT.TFile(indir+'/'+sd+'/scan_'+sd+'.root', 'read')
            tree = f.Get('fitresults')
     
            vals = []
            for ev in tree:
                vals.append( [getattr(ev, par), 2.*ev.nllval  ] )
            vals = sorted(vals)
            lvals = vals[:len(vals)/2]
            rvals = vals[len(vals)/2:]
     
            graph = ROOT.TGraph(len(vals), array.array('d', [x[1] for x in vals]), array.array('d', [y[0] for y in vals]) )
     
            best = graph.Eval(0.)
            lgraph = ROOT.TGraph(len(lvals), array.array('d', [x[1] for x in lvals]), array.array('d', [y[0] for y in lvals]) )
            rgraph = ROOT.TGraph(len(rvals), array.array('d', [x[1] for x in rvals]), array.array('d', [y[0] for y in rvals]) )
            sol1  = lgraph.Eval(1.)
            sol2  = rgraph.Eval(1.)
     
            _dict[par] = (best, sol1, sol2)
     
        return _dict


    def effSigma(self, histo):
        xaxis = histo.GetXaxis()
        nb = xaxis.GetNbins()
        xmin = xaxis.GetXmin()
        ave = histo.GetMean()
        rms = histo.GetRMS()
        total=histo.Integral()
        if total < 100: 
            print("effsigma: Too few entries to compute it: ", total)
            return 0.
        ierr=0
        ismin=999
        rlim=0.683*total
        bwid = xaxis.GetBinWidth(1)
        nrms=int(rms/bwid)
        if nrms > nb/10: nrms=int(nb/10) # Could be tuned...
        widmin=9999999.
        for iscan in xrange(-nrms,nrms+1): # // Scan window centre 
            ibm=int((ave-xmin)/bwid)+1+iscan
            x=(ibm-0.5)*bwid+xmin
            xj=x; xk=x;
            jbm=ibm; kbm=ibm;
            bin=histo.GetBinContent(ibm)
            total=bin
            for j in xrange(1,nb):
                if jbm < nb:
                    jbm += 1
                    xj += bwid
                    bin=histo.GetBinContent(jbm)
                    total += bin
                    if total > rlim: break
                else: ierr=1
                if kbm > 0:
                    kbm -= 1
                    xk -= bwid
                    bin=histo.GetBinContent(kbm)
                    total+=bin
                if total > rlim: break
                else: ierr=1
            dxf=(total-rlim)*bwid/bin
            wid=(xj-xk+bwid-dxf)*0.5
            if wid < widmin:
                widmin=wid
                ismin=iscan
        if ismin == nrms or ismin == -nrms: ierr=3
        if ierr != 0: print("effsigma: Error of type ", ierr)
        return widmin

    def doShadedUncertainty(self,h):
        xaxis = h.GetXaxis()
        points = []; errors = []
        for i in xrange(h.GetNbinsX()):
            N = h.GetBinContent(i+1); dN = h.GetBinError(i+1);
            if N == 0 and dN == 0: continue
            x = xaxis.GetBinCenter(i+1);
            points.append( (x,N) )
            EYlow, EYhigh  = dN, min(dN,N);
            EXhigh, EXlow = (xaxis.GetBinUpEdge(i+1)-x, x-xaxis.GetBinLowEdge(i+1))
            errors.append( (EXlow,EXhigh,EYlow,EYhigh) )
        ret = ROOT.TGraphAsymmErrors(len(points))
        ret.SetName(h.GetName()+"_errors")
        for i,((x,y),(EXlow,EXhigh,EYlow,EYhigh)) in enumerate(zip(points,errors)):
            ret.SetPoint(i, x, y)
            ret.SetPointError(i, EXlow,EXhigh,EYlow,EYhigh)
        ret.SetFillStyle(3244);
        ret.SetFillColor(ROOT.kGray+2)
        ret.SetMarkerStyle(0)
        ret.Draw("PE2 SAME")
        return ret

    def safecolor(self, index):
        SAFE_COLOR_LIST=[ROOT.kBlack, ROOT.kRed, ROOT.kGreen+2, ROOT.kBlue, ROOT.kMagenta+1, ROOT.kOrange+7, ROOT.kCyan+1, ROOT.kGray+2, ROOT.kViolet+5, ROOT.kSpring+5, ROOT.kAzure+1, ROOT.kPink+7, ROOT.kOrange+3, ROOT.kBlue+3, ROOT.kMagenta+3, ROOT.kRed+2]+range(11,40)
        if index<len(SAFE_COLOR_LIST): return SAFE_COLOR_LIST[index]
        else: return index

    def getCoeffs(self,xL,xR,x0,err_xL,err_xR,err_x0, toyEvents=10000):        
        histos = { 'a0': ROOT.TH1D('a0','',100,-0.4,0.6),
                   'a4': ROOT.TH1D('a4','',100,-0.4,3.0)
                   }
        # given that we use mean and RMS, for security reason use Under/Overflow values to compute them
        # so we don't have to tune the range (even though the one about is already very sensible for a0 and a4)
        histos['a0'].SetStatOverflows(1)
        histos['a4'].SetStatOverflows(1)
        #print("getCoeffs: ",toyEvents," toyMC running...")
        for i in xrange(toyEvents):
            ixL = np.random.normal(xL,err_xL)
            ixR = np.random.normal(xR,err_xR)
            ix0 = np.random.normal(x0,err_x0)
            sumPol = ixL+ixR+ix0
            histos['a0'].Fill(2*ix0/sumPol)
            #histos['a0'].Fill(2*(1-ixL-ixR)/sumPol)
            histos['a4'].Fill(2*(ixL-ixR)/sumPol)
        #print("toyMC done")
        ret = {}
        for k,h in histos.iteritems():
            ret[k] = (h.GetMean(),h.GetRMS())
        return ret
            
    def getChargeAsy(self,xplus,xminus,err_xplus,err_xminus,toyEvents=10000):
        histo = ROOT.TH1F('hasy','',100,-0.05,0.5)
        for i in xrange(toyEvents):
            ixplus  = np.random.normal(xplus,err_xplus)
            ixminus = np.random.normal(xminus,err_xminus)
            histo.Fill((ixplus-ixminus)/(ixplus+ixminus))
        #print("toyMC done")
        ret = {'asy': (histo.GetMean(),histo.GetRMS())}
        del histo
        return ret

    def getChargeAsyFromTH1pair(self,h1, h2, toyEvents=10000, name='asy'):
        # assuming h1 and h2 have same binning
        print("Inside getChargeAsyFromTH1pair()")
        binsx = [h1.GetXaxis().GetBinLowEdge(i) for i in range(1,2+h1.GetNbinsX())]
        histo = ROOT.TH1D(name,'',len(binsx)-1, array.array('d',binsx))
        histotmp = ROOT.TH1D(name+'tmp','',200,-1.0,1.0)
        for i in range(1,1+h1.GetNbinsX()):
            histotmp.Reset("ICESM")
            for j in xrange(toyEvents):
                ixplus  = np.random.normal(h1.GetBinContent(i),h1.GetBinError(i))
                ixminus = np.random.normal(h2.GetBinContent(i),h2.GetBinError(i))
                histotmp.Fill((ixplus-ixminus)/(ixplus+ixminus))
            histo.SetBinContent(i, histotmp.GetMean()  )
            histo.SetBinError(  i, histotmp.GetStdDev())
        #print("toyMC done")
        return histo


    def getNEffStat(self, s):
        a = s.split('EffStat')[1]
        a = a.replace('minus','').replace('plus','')
        a = a.replace('mu','').replace('el','')
        return int(a)

    def getNFromString(self, s, chooseIndex=0, useAll=False):
        los = [ int(i) for i in re.findall(r'\d+', s) ]
        if len(los) == 0: return 0
        if len(los) == 1: return los[0]
        if len(los) > 1:
            if useAll:                
                return tuple(los)
            else:
                return los[min(chooseIndex,len(los)-1)]
        return 0
        
    def getChannelFromFitresults(self, fitresults):
        t_fitres = fitresults.Get('fitresults')
        channel = 'el' if 'ErfPar0EffStat1elminus' in t_fitres.GetListOfBranches() else 'mu'
        return channel

    def wxsec(self,generator='mcatnlo'):
        # FEWZ3.1 comes from https://twiki.cern.ch/twiki/bin/viewauth/CMS/StandardModelCrossSectionsat13TeV
        # 770.9  pb is the PDF uncertainty
        # +165.7 -88.2 is the QCD scales uncertainty
        # for the native mcatnlo 110 pb is added in quadrature since it is the stat error on the partial sample (JB recomputed it -cf mail of 8/8/2019)
        xsec = {'fewz3p1': (3*20508.9,math.hypot(165.7,770.9),math.hypot(88.2,770.9)),
                'mcatnlo': (60400,    math.hypot(110,math.hypot(165.7,770.9)),math.hypot(110,math.hypot(88.2,770.9)))}
        if generator not in xsec:
            print("ERROR! Generator ",generator," unknown. Returning 0 xsec")
            return (0,-1,-1)
        return xsec[generator]

    def getL1SF(self,pt,eta,histo):
        ret = (0,0)
        if not histo: 
            print("The ", histo, " is not present in the file ",rfile)
            return ret
        etabin = max(1, min(histo.GetNbinsX(), histo.GetXaxis().FindFixBin(eta)))
        ptbin  = max(1, min(histo.GetNbinsY(), histo.GetYaxis().FindFixBin(pt)))
        ret = ( histo.GetBinContent(etabin,ptbin), histo.GetBinError(etabin,ptbin) )
        return ret

    def getEffSyst(self,pdgId):
        if abs(pdgId)==11:
            etabins_el = array.array('f',[0.0, 1.0, 1.479, 2.0, 2.2, 2.5])
            h1d = ROOT.TH1F('effsysth','',len(etabins_el)-1,etabins_el)
            h1d.SetBinContent(1, 0.006)
            h1d.SetBinContent(2, 0.008)
            h1d.SetBinContent(3, 0.013)
            h1d.SetBinContent(4, 0.016)
            h1d.SetBinContent(5, 0.019)
        elif abs(pdgId)==13:
            etabins_mu = array.array('f',[0.0, 1.0, 1.5, 2.4])
            h1d = ROOT.TH1F('effsysth','',len(etabins_mu)-1,etabins_mu)
            h1d.SetBinContent(1, 0.002)
            h1d.SetBinContent(2, 0.004)
            h1d.SetBinContent(3, 0.014)
        return h1d

    def getExclusiveBinnedSyst(self,th1):
        th1_excl = th1.Clone(th1.GetName()+'_exclusive')
        for ibin in range(1,th1.GetNbinsX()+1):
            tmp_val_sq = math.pow(th1.GetBinContent(ibin),2)
            for inner_bin in range(1,ibin):
                tmp_val_sq = tmp_val_sq - math.pow(th1_excl.GetBinContent(inner_bin),2)
                if tmp_val_sq<0:
                    print("WARNING! getExclusiveBinnedSyst cropping syst at 0. ")
                    tmp_val_sq = 0
                    break
            th1_excl.SetBinContent(ibin,math.sqrt(tmp_val_sq))
        return th1_excl

    def getRochesterUncertainty(self,charge,etabin,ptbin,syst_histo,averagept):
        ## the syst 3 is the only one with a pt-shape, so maintain it
        return np.mean(syst_histo, axis=1)[etabin] if averagept else syst_histo[etabin,ptbin]

