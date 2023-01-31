
import sys,copy,array,os,subprocess,math
import ROOT


if not 'CMSSW_BASE' in os.environ:
    print("CMSSW/combine not found")
    quit()
    
cmssw_base = os.environ['CMSSW_BASE']
combine_base = "{cmssw_base}/src/HiggsAnalysis/CombinedLimit/scripts/".format(cmssw_base=cmssw_base)

class Combine():

    def __init__(self, cardName, outDir):
    
        self.card = os.path.realpath(cardName)
        self.runDir = os.path.dirname(os.path.realpath(card))
        self.outDir = outDir
        
        self.hdf5 = self.card.replace(".txt", ".hdf5")
        self.out = self.card.replace(".txt", "_output.root")

    def text2hdf5(self):
    
        print("### Run text2hdf5")
        cmd = "python {combine_base}/text2hdf5.py {card} -o {hdf5} --X-allow-no-signal --X-allow-no-background".format(combine_base=combine_base, card=self.card, hdf5=self.hdf5)
        subprocess.call(cmd, shell=True, cwd=self.runDir)

    def run_differential(self):

        print("### Run differential fit")
        cmd = "python {combine_base}/combinetf.py {hdf5} -o {out} -t -1 --binByBinStat --expectSignal=1 --doImpacts --saveHists --computeHistErrors".format(combine_base=combine_base, out=self.out, hdf5=self.hdf5)
        subprocess.call(cmd, shell=True, cwd=self.runDir)
    
    
def run_differential_norm():

    pass


if __name__ == "__main__":

    card = "CombineStudies/lowPU_differential/lowPU_Zmumu_RawPFMET_differential.txt"
    outDir = "/eos/user/j/jaeyserm/www/wmass/combine/lowPU_differential/"
    c = Combine(card, outDir)
    c.text2hdf5()
    c.run_differential()
    
    
    