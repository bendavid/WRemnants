
import sys,copy,array,os,subprocess,math
import ROOT


if not 'CMSSW_BASE' in os.environ:
    print("CMSSW/combine not found")
    quit()
    
cmssw_base = os.environ['CMSSW_BASE']
combine_base = "{cmssw_base}/src/HiggsAnalysis/CombinedLimit/scripts/".format(cmssw_base=cmssw_base)

class Combine():

    def __init__(self, cardName, outDir, norm=False):
    
        self.norm = norm
        self.card = os.path.realpath(cardName)
        self.runDir = os.path.dirname(os.path.realpath(card))
        self.outDir = outDir
        
        self.hdf5 = self.card.replace(".txt", ".hdf5")
        self.out = self.card.replace(".txt", "_output.root")
        self.cat = "mumu"
        
        if norm:
        
            self.card_bare = self.card
            self.card_xsec = self.card.replace(".txt", "_xsec.txt")
            self.card = self.card.replace(".txt", "_tot.txt")
            self.hdf5 = self.card.replace(".txt", ".hdf5")
            self.out = self.card.replace(".txt", "_output.root")
        
            print("### Run combineCards")
            cmd = "python {combine_base}/combineCards.py --noDirPrefix {cat}={card_bare} {cat}_xsec={card_xsec} > {card}".format(combine_base=combine_base, card=self.card, card_bare=self.card_bare, card_xsec=self.card_xsec, cat=self.cat)
            subprocess.call(cmd, shell=True, cwd=self.runDir)
            


    def text2hdf5(self):
    
        print("### Run text2hdf5")
        opts = ""
        if self.norm:
            opts = "--maskedChan={cat}_xsec".format(cat=self.cat)
            
        cmd = "python {combine_base}/text2hdf5.py {card} -o {hdf5} --X-allow-no-signal --X-allow-no-background {opts}".format(combine_base=combine_base, card=self.card, hdf5=self.hdf5, opts=opts)
        subprocess.call(cmd, shell=True, cwd=self.runDir)

    def run_differential(self):

        self.text2hdf5()

        print("### Run differential fit")
        opts = ""

        cmd = "python {combine_base}/combinetf.py {hdf5} -o {out} -t -1 --binByBinStat --expectSignal=1 --doImpacts --saveHists --computeHistErrors {opts}".format(combine_base=combine_base, out=self.out, hdf5=self.hdf5, opts=opts)
        subprocess.call(cmd, shell=True, cwd=self.runDir)
        
    



if __name__ == "__main__":

    card = "CombineStudies/lowPU_differential/lowPU_Zmumu_RawPFMET_differential_stat.txt"
    #card = "CombineStudies/lowPU_differential/lowPU_Zmumu_RawPFMET_differential.txt"
    outDir = "/eos/user/j/jaeyserm/www/wmass/combine/lowPU_differential/"
    c = Combine(card, outDir, norm=True)
    c.run_differential()
    
    
    