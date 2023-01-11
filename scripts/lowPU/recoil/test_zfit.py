#  Copyright (c) 2022 zfit

import numpy as np
import time

import sys
sys.path.insert(0, "/home/j/jaeyserm/analysis/WRemnants/env/lib/python3.10/site-packages")
import zfit

#t = zfit.run.n_cpu()
#print(t)

zfit.run.set_n_cpu(1)
zfit.run.set_cpus_explicit(1, 1)

zfit.run.set_n_cpu(256)
zfit.run.set_cpus_explicit(256, 256)
print(zfit.run.n_cpu)

#zfit.run.set_autograd_mode(True)
#zfit.run.set_graph_mode(True)


# create space
obs = zfit.Space("x", limits=(-50, 50))
obs_binned = obs.with_binning(100)

# parameters
mu1 = zfit.Parameter("mu1", 10, -50, 50)
sigma1 = zfit.Parameter("sigma1", 5, 0.5, 10)
gauss1 = zfit.pdf.Gauss(mu=mu1, sigma=sigma1, obs=obs, name="gauss1")
gauss1_binned = zfit.pdf.BinnedFromUnbinnedPDF(gauss1, obs_binned)

mu2 = zfit.Parameter("mu2", 0, -50, 50)
sigma2 = zfit.Parameter("sigma2", 3, 0.5, 10)
gauss2 = zfit.pdf.Gauss(mu=mu2, sigma=sigma2, obs=obs, name="gauss2")
gauss2_binned = zfit.pdf.BinnedFromUnbinnedPDF(gauss2, obs_binned)

mu3 = zfit.Parameter("mu3", -10, -50, 50)
sigma3 = zfit.Parameter("sigma3", 5, 0.5, 10)
gauss3 = zfit.pdf.Gauss(mu=mu3, sigma=sigma3, obs=obs, name="gauss3")
gauss3_binned = zfit.pdf.BinnedFromUnbinnedPDF(gauss3, obs_binned)

mu4 = zfit.Parameter("mu4", 30, -50, 50)
sigma4 = zfit.Parameter("sigma4", 1, 0.5, 10)
gauss4 = zfit.pdf.Gauss(mu=mu4, sigma=sigma4, obs=obs, name="gauss4")
gauss4_binned = zfit.pdf.BinnedFromUnbinnedPDF(gauss4, obs_binned)

f1 = zfit.Parameter("f1", 0.2, 0, 1)
f2 = zfit.Parameter("f2", 0.5, 0, 1)
f3 = zfit.Parameter("f3", 0.2, 0, 1)

f1.floating = True
f2.floating = True
f3.floating = True

gauss = zfit.pdf.SumPDF([gauss1, gauss2, gauss3, gauss4], fracs=[f1, f2, f3], obs=obs, name="gauss")

minimizer = zfit.minimize.Minuit(mode=2)
#minimizer = zfit.minimize.ScipyTrustConstrV1() 

binned = True
if binned:

    # model building, pdf creation
    #gauss_binned = gauss1_binned
    #gauss_binned = zfit.pdf.BinnedFromUnbinnedPDF(gauss, obs_binned) # 10 seconds
    gauss_binned = zfit.pdf.BinnedSumPDF([gauss1_binned, gauss2_binned, gauss3_binned, gauss4_binned], fracs=[f1, f2, f3], obs=obs_binned, name="gauss") # 5 seconds

    data_binned = gauss_binned.sample(n=1000000)

    # create NLL
    nll = zfit.loss.BinnedNLL(model=gauss_binned, data=data_binned)
    print("Start miminization")
    start = time.time()
    result = minimizer.minimize(nll)
    end = time.time()
    print("##### MIMIMIZATION TIME", (end-start))
    print(result)


    
else:

    data = gauss.create_sampler(n=300000)
    data.resample()
    nll = zfit.loss.UnbinnedNLL(model=gauss, data=data)
    print("Start miminization")
    start = time.time()
    result = minimizer.minimize(nll)
    end = time.time()
    print("##### MIMIMIZATION TIME", (end-start))
    print(result)