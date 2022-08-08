import lz4.frame
import pickle

def read_and_scale(fname, proc, histname):
    with lz4.frame.open(fname) as f:
        out = pickle.load(f)
        
    return out[proc]["output"][histname]*out[proc]["dataset"]["xsec"]/out[proc]["weight_sum"]
