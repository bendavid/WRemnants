
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
        self.runDir = os.path.dirname(os.path.realpath(self.card))
        self.outDir = outDir
        
        self.hdf5 = self.card.replace(".txt", ".hdf5")
        self.out = self.card.replace(".txt", "_output.root")
        self.cat = "mumu"
        


    def text2hdf5(self):
    
        print("### Run text2hdf5")
        opts = ""

        cmd = "python {combine_base}/text2hdf5.py {card} -o {hdf5} --X-allow-no-signal --X-allow-no-background {opts}".format(combine_base=combine_base, card=self.card, hdf5=self.hdf5, opts=opts)
        subprocess.call(cmd, shell=True, cwd=self.runDir)

    def run_fit(self):

        self.text2hdf5()

        print("### Run fit")
        opts = ""

        cmd = "python {combine_base}/combinetf.py {hdf5} -o {out} -t -1 --binByBinStat --expectSignal=1 --doImpacts --saveHists --computeHistErrors {opts}".format(combine_base=combine_base, out=self.out, hdf5=self.hdf5, opts=opts)
        subprocess.call(cmd, shell=True, cwd=self.runDir)
    

def combineCards(card1, card2, cardOut):

    card1 = os.path.realpath(card1)
    card2 = os.path.realpath(card2)
    subprocess.call("touch {cardOut}".format(cardOut=cardOut), shell=True)
    cardOut = os.path.realpath(cardOut)
    runDir = os.path.dirname(os.path.realpath(cardOut))

    cmd = "python {combine_base}/combineCards.py --noDirPrefix {card1} {card2} > {cardOut}".format(combine_base=combine_base, card1=card1, card2=card2, cardOut=cardOut)
    subprocess.call(cmd, shell=True, cwd=runDir)
    



if __name__ == "__main__":

    outDir = "/eos/user/j/jaeyserm/www/wmass/combine/lowPU_highPU/"
    
    card_plus = "ZMassWLike/ZMassWLike_plus.txt"
    
    c = Combine(card_plus, outDir)
    c.run_fit()
    
    card_minus = "ZMassWLike/ZMassWLike_minus.txt"
    c = Combine(card_minus, outDir)
    c.run_fit()
    
    card_z = "ZMassWLike/ZMassWLike.txt"
    combineCards(card_plus, card_minus, card_z)
    c = Combine(card_z, outDir)
    c.run_fit()
    
    card_lowPU_wmass = "CombineStudies/lowPU_wmass/lowPU_Zmumu_RawPFMET_wmass.txt"
    card_combined = "CombineStudies/lowPU_wmass/lowPU_highPU.txt"
    combineCards(card_z, card_lowPU_wmass, card_combined)
    c = Combine(card_combined, outDir)
    c.run_fit()
    
