

import math
import tensorflow as tf
import baseFunctions as bf


    

def ewk_para_qT(xvals, p, *args):

    mean1 = p[3]
    mean2 = p[4]
    mean3 = p[5]
    
    gauss1 = bf.func_gauss(xvals[0], mean1, p[0])
    gauss2 = bf.func_gauss(xvals[0], mean2, p[1])
    gauss3 = bf.func_gauss(xvals[0], mean3, p[2])

    n1 = 0.33
    n2 = 0.33
    n3 = tf.constant(1.0, dtype=tf.float64) - n2 - n1

    pdf = n1*gauss1 + n2*gauss2 + n3*gauss3
    return pdf
    
def ewk_para_qT_cond(xvals, p, *args):
    '''
    sigma1 = quadratic(xvals[0], p[0], p[1], p[2])
    sigma2 = quadratic(xvals[0], p[3], p[4], p[5])
    sigma3 = quadratic(xvals[0], p[6], p[7], p[8])
    
    #pp = tf.constant([2.99723e+00, -1.47337e-01, 1.61576e-02, -1.42718e-04, -9.91942e-07, 1.47306e-08, -3.50397e-11], dtype=tf.float64)
    #mean = pol6(xvals[0], pp[0], pp[1], pp[2], pp[3], pp[4], pp[5], pp[6])
    mean1 = linear(xvals[0], p[9], p[10])
    mean2 = linear(xvals[0], p[11], p[2])
    mean3 = linear(xvals[0], p[13], p[14])

    gauss1 = func_gauss(xvals[1], mean1, sigma1)
    gauss2 = func_gauss(xvals[1], mean2, sigma2)
    gauss3 = func_gauss(xvals[1], mean3, sigma3)

    n1 = p[15] 
    n2 = p[16]
    n3 = tf.constant(1, dtype=tf.float64) - n2 - n1
    '''
    sigma1 = bf.power(xvals[0], p[0], p[1], p[2])
    sigma2 = bf.power(xvals[0], p[3], p[4], p[5])
    sigma3 = bf.power(xvals[0], p[6], p[7], p[8])
    #sigma1 = quadratic(xvals[0], p[0], p[1], p[2])
    #sigma2 = quadratic(xvals[0], p[3], p[4], p[5])
    #sigma3 = quadratic(xvals[0], p[6], p[7], p[8])
    
    #pp = tf.constant([2.19057e+00, -2.37865e-01, 4.02425e-02, -1.34963e-03, 2.19717e-05, -1.73520e-07, 5.31624e-10], dtype=tf.float64)
    #mean = pol6(xvals[0], pp[0], pp[1], pp[2], pp[3], pp[4], pp[5], pp[6])
    #mean1 = quadratic(xvals[0], p[9], p[10], p[11])
    #mean2 = quadratic(xvals[0], p[12], p[13], p[14])
    #mean3 = quadratic(xvals[0], p[15], p[16], p[17])
    mean1 = bf.linear(xvals[0], p[9], p[10])
    mean2 = bf.linear(xvals[0], p[11], p[12])
    mean3 = bf.linear(xvals[0], p[13], p[14])

    gauss1 = bf.func_gauss(xvals[1], mean1, sigma1)
    gauss2 = bf.func_gauss(xvals[1], mean2, sigma2)
    gauss3 = bf.func_gauss(xvals[1], mean3, sigma3)

    #n1 = p[12]
    #n2 = p[13]
    #n3 = tf.constant(1, dtype=tf.float64) - n2 - n1
    n1 = p[15]
    n2 = p[16]
    n3 = 1.0 - n2 - n1
    
    pdf = n1*gauss1 + n2*gauss2 + n3*gauss3
    return pdf
 

 
def ewk_perp(xvals, p, *args):
    
    gauss1 = bf.func_gauss(xvals[0], tf.constant(0, dtype=tf.float64), p[0])
    gauss2 = bf.func_gauss(xvals[0], tf.constant(0, dtype=tf.float64), p[1])
    gauss3 = bf.func_gauss(xvals[0], tf.constant(0, dtype=tf.float64), p[2])

    n1 = 0.76
    n2 = 0.22
    n3 = tf.constant(1.0, dtype=tf.float64) - n2 - n1

    pdf = n1*gauss1 + n2*gauss2 + n3*gauss3
    return pdf
 
def ewk_perp_cond(xvals, p, *args):

    sigma1 = bf.power(xvals[0], p[0], p[1], p[2])
    sigma2 = bf.power(xvals[0], p[3], p[4], p[5])
    sigma3 = bf.linear(xvals[0], p[6], p[7])

    gauss1 = bf.func_gauss(xvals[1], tf.constant(0, dtype=tf.float64), sigma1)
    gauss2 = bf.func_gauss(xvals[1], tf.constant(0, dtype=tf.float64), sigma2)
    gauss3 = bf.func_gauss(xvals[1], tf.constant(0, dtype=tf.float64), sigma3)

    n1 = 0.76
    n2 = 0.22
    n3 = tf.constant(1, dtype=tf.float64) - n2 - n1

    pdf = n1*gauss1 + n2*gauss2 + n3*gauss3
    return pdf
 
 
 
 
def ttbar_para_qT(xvals, p, *args):

    gauss1 = bf.func_gauss(xvals[0], p[3], p[0])
    gauss2 = bf.func_gauss(xvals[0], p[4], p[1])
    gauss3 = bf.func_gauss(xvals[0], p[5], p[2])

    n1 = tf.constant(0.65, dtype=tf.float64)
    n2 = tf.constant(0.25, dtype=tf.float64)
    n3 = tf.constant(1.0, dtype=tf.float64) - n2 - n1

    pdf = n1*gauss1 + n2*gauss2 + n3*gauss3
    return pdf
     
def ttbar_para_qT_cond(xvals, p, *args):

    sigma1 = bf.quadratic(xvals[0], p[0], p[1], p[2])
    sigma2 = bf.quadratic(xvals[0], p[3], p[4], p[5])
    sigma3 = bf.quadratic(xvals[0], p[6], p[7], p[8])
    
    mean1 = bf.quadratic(xvals[0], p[9], p[10], p[11])
    mean2 = bf.quadratic(xvals[0], p[12], p[13], p[14])
    mean3 = bf.quadratic(xvals[0], p[15], p[16], p[17])

    gauss1 = bf.func_gauss(xvals[1], mean1, sigma1)
    gauss2 = bf.func_gauss(xvals[1], mean2, sigma2)
    gauss3 = bf.func_gauss(xvals[1], mean3, sigma3)

    n1 = p[18] 
    n2 = p[19]
    n3 = tf.constant(1, dtype=tf.float64) - n2 - n1

    pdf = n1*gauss1 + n2*gauss2 + n3*gauss3
    return pdf
 
 
 
def ttbar_perp(xvals, p, *args):
    
    gauss1 = bf.func_gauss(xvals[0], tf.constant(0, dtype=tf.float64), p[0])
    gauss2 = bf.func_gauss(xvals[0], tf.constant(0, dtype=tf.float64), p[1])
    gauss3 = bf.func_gauss(xvals[0], tf.constant(0, dtype=tf.float64), p[2])

    n1 = tf.constant(0.5, dtype=tf.float64)
    n2 = tf.constant(0.4, dtype=tf.float64)
    n3 = tf.constant(1.0, dtype=tf.float64) - n2 - n1

    pdf = n1*gauss1 + n2*gauss2 + n3*gauss3
    return pdf
    
def ttbar_perp_cond(xvals, p, *args):

    sigma1 = bf.quadratic(xvals[0], p[0], p[1], p[2])
    sigma2 = bf.quadratic(xvals[0], p[3], p[4], p[5])
    sigma3 = bf.quadratic(xvals[0], p[6], p[7], p[8])

    gauss1 = bf.func_gauss(xvals[1], tf.constant(0, dtype=tf.float64), sigma1)
    gauss2 = bf.func_gauss(xvals[1], tf.constant(0, dtype=tf.float64), sigma2)
    gauss3 = bf.func_gauss(xvals[1], tf.constant(0, dtype=tf.float64), sigma3)

    n1 = p[9] #tf.constant(0.5, dtype=tf.float64)
    n2 = p[10] #tf.constant(0.49, dtype=tf.float64)
    n3 = tf.constant(1, dtype=tf.float64) - n2 - n1

    pdf = n1*gauss1 + n2*gauss2 + n3*gauss3
    return pdf
 



def dy_para_qT(xvals, p, *args):

    mean1 = p[4]
    mean2 = p[5]
    mean3 = p[6]
    #mean4 = p[7]
    
    gauss1 = bf.func_gauss(xvals[0], mean1, p[0])
    gauss2 = bf.func_gauss(xvals[0], mean2, p[1])
    gauss3 = bf.func_gauss(xvals[0], mean3, p[2])
    gauss4 = bf.func_gauss(xvals[0], mean3, p[3])

    #n1 = tf.constant(0.01, dtype=tf.float64)
    #n2 = tf.constant(0.2, dtype=tf.float64)
    #n3 = tf.constant(0.6, dtype=tf.float64)
    #n4 = tf.constant(1.0, dtype=tf.float64) - n1 - n2 - n3
    
    n1 = tf.constant(0.1, dtype=tf.float64)
    n2 = tf.constant(0.2, dtype=tf.float64)
    n3 = tf.constant(0.25, dtype=tf.float64)
    n4 = tf.constant(1.0, dtype=tf.float64) - n1 - n2 - n3

    pdf = n1*gauss1 + n2*gauss2 + n3*gauss3 + n4*gauss4
    return pdf

def dy_para_qT_cond(xvals, p, *args):

    sigma1 = bf.power(xvals[0], p[0], p[1], p[2])
    sigma2 = bf.power(xvals[0], p[3], p[4], p[5])
    sigma3 = bf.power(xvals[0], p[6], p[7], p[8])
    sigma4 = bf.power(xvals[0], p[9], p[10], p[11])
    
    mean1 = bf.power(xvals[0], p[12], p[13], p[14])
    mean2 = bf.power(xvals[0], p[15], p[16], p[17])
    mean3 = bf.power(xvals[0], p[18], p[19], p[20])
    #mean4 = bf.power(xvals[0], p[21], p[22], p[23])
    
    gauss1 = bf.func_gauss(xvals[1], mean1, sigma1)
    gauss2 = bf.func_gauss(xvals[1], mean2, sigma2)
    gauss3 = bf.func_gauss(xvals[1], mean3, sigma3)
    gauss4 = bf.func_gauss(xvals[1], mean3, sigma4)

    n1 = p[21]
    n2 = p[22]
    n3 = p[23]
    n4 = tf.constant(1.0, dtype=tf.float64) - n1 - n2 - n3
   
    pdf = n1*gauss1 + n2*gauss2 + n3*gauss3 + n4*gauss4
    return pdf
 
 
 
def dy_perp(xvals, p, *args):
   
    gauss1 = bf.func_gauss(xvals[0], p[4], p[0])
    gauss2 = bf.func_gauss(xvals[0], p[5], p[1])
    gauss3 = bf.func_gauss(xvals[0], p[4], p[2])
    gauss4 = bf.func_gauss(xvals[0], p[5], p[3])
    
    n1 = tf.constant(0.1, dtype=tf.float64)
    n2 = tf.constant(0.2, dtype=tf.float64)
    n3 = tf.constant(0.25, dtype=tf.float64)
    n4 = tf.constant(1.0, dtype=tf.float64) - n1 - n2 - n3

    pdf = n1*gauss1 + n2*gauss2 + n3*gauss3 + n4*gauss4
    return pdf
    
def dy_perp_cond(xvals, p, *args):

    sigma1 = bf.power(xvals[0], p[0], p[1], p[2])
    sigma2 = bf.power(xvals[0], p[3], p[4], p[5])
    sigma3 = bf.power(xvals[0], p[6], p[7], p[8])
    sigma4 = bf.power(xvals[0], p[9], p[10], p[11])
    
    mean1 = bf.linear(xvals[0], p[12], p[13])
    mean2 = bf.linear(xvals[0], p[14], p[15])
    
    gauss1 = bf.func_gauss(xvals[1], mean1, sigma1)
    gauss2 = bf.func_gauss(xvals[1], mean2, sigma2)
    gauss3 = bf.func_gauss(xvals[1], mean1, sigma3)
    gauss4 = bf.func_gauss(xvals[1], mean2, sigma4)

    n1 = p[16]
    n2 = p[17]
    n3 = p[18]
    n4 = tf.constant(1.0, dtype=tf.float64) - n1 - n2 - n3
   
    pdf = n1*gauss1 + n2*gauss2 + n3*gauss3 + n4*gauss4
    return pdf
    

 

 
def data_para_qT(xvals, p, *args):

    gauss1 = bf.func_gauss(xvals[0], p[3], p[0])
    gauss2 = bf.func_gauss(xvals[0], p[4], p[1])
    gauss3 = bf.func_gauss(xvals[0], p[4], p[2])

    n1 = tf.constant(0.2, dtype=tf.float64)
    n2 = tf.constant(0.12, dtype=tf.float64)
    n3 = tf.constant(1.0, dtype=tf.float64) - n1 - n2
   
    pdf_sig = n1*gauss1 + n2*gauss2 + n3*gauss3
    
    pdf_ttbar = ttbar_para_qT_cond([args[0]['qT'], xvals[0]], args[0]['parms_ttbar'])
    pdf_ewk = ewk_para_qT_cond([args[0]['qT'], xvals[0]], args[0]['parms_ewk'])
    
    qTbin = tf.cast(tf.math.multiply(args[0]['qT'], tf.constant(2, dtype=tf.float64)), dtype=tf.int64)
    k1 = tf.gather(args[0]['fracs_ttbar'], qTbin)
    k2 = tf.gather(args[0]['fracs_ewk'], qTbin)
    k3 = tf.constant(1.0, dtype=tf.float64) - k1 - k2
    pdf = k1*pdf_ttbar + k2*pdf_ewk + k3*pdf_sig
    return pdf
 
def data_para_qT_cond_mctemplates(xvals, p, *args):

    pdf_ttbar = ttbar_para_qT_cond(xvals, args[0]['parms_ttbar'])
    pdf_ewk = ewk_para_qT_cond(xvals, args[0]['parms_ewk'])
    
    
    sigma1 = bf.power(xvals[0], p[0], p[1], p[2])
    sigma2 = bf.power(xvals[0], p[3], p[4], p[5])
    sigma3 = bf.power(xvals[0], p[6], p[7], p[8])
    sigma4 = bf.power(xvals[0], p[9], p[10], p[11])
    
    mean1 = bf.power(xvals[0], p[12], p[13], p[14])
    mean2 = bf.power(xvals[0], p[15], p[16], p[17])
    #mean3 = bf.power(xvals[0], p[18], p[19], p[20])
    #mean4 = power(xvals[0], p[21], p[22], p[23])
    
    gauss1 = bf.func_gauss(xvals[1], mean1, sigma1)
    gauss2 = bf.func_gauss(xvals[1], mean1, sigma2)
    gauss3 = bf.func_gauss(xvals[1], mean2, sigma3)
    gauss4 = bf.func_gauss(xvals[1], mean2, sigma4)

    n1 = p[18]
    n2 = p[19]
    n3 = 0.5886822304617599 # p[22]
    n4 = tf.constant(1.0, dtype=tf.float64) - n1 - n2 - n3
   
    pdf_sig = n1*gauss1 + n2*gauss2 + n3*gauss3 + n4*gauss4
    
    # LnN: https://www.physics.ucla.edu/~cousins/stats/cousins_lognormal_prior.pdf
    # case needed as p[x] is double during plotting?!
    #ttbar_unc = math.log(1.001)
    #ttbar_constr = func_lognormal(tf.cast(p[11], dtype=tf.float64), tf.constant(0.0, dtype=tf.float64), tf.constant(ttbar_unc, dtype=tf.float64))
    
    #ewk_unc = math.log(1.01)
    #ewk_constr = func_lognormal(tf.cast(p[12], dtype=tf.float64), tf.constant(0.0, dtype=tf.float64), tf.constant(ewk_unc, dtype=tf.float64))

    qTbin = tf.cast(tf.math.multiply(xvals[0], tf.constant(2, dtype=tf.float64)), dtype=tf.int64)
    k1 = tf.gather(args[0]['fracs_ttbar'], qTbin)
    k2 = tf.gather(args[0]['fracs_ewk'], qTbin)
    k3 = tf.constant(1.0, dtype=tf.float64) - k1 - k2
    pdf = k1*pdf_ttbar + k2*pdf_ewk + k3*pdf_sig
    return pdf
    
 
def data_para_qT_cond(xvals, p, *args):
 
    sigma1 = bf.power(xvals[0], p[0], p[1], p[2])
    sigma2 = bf.power(xvals[0], p[3], p[4], p[5])
    sigma3 = bf.power(xvals[0], p[6], p[7], p[8])
    
    mean1 = bf.power(xvals[0], p[9], p[10], p[11])
    mean2 = bf.power(xvals[0], p[12], p[13], p[14])
    
    gauss1 = bf.func_gauss(xvals[1], mean1, sigma1)
    gauss2 = bf.func_gauss(xvals[1], mean2, sigma2)
    gauss3 = bf.func_gauss(xvals[1], mean2, sigma3)

    n1 = p[15]
    n2 = p[16]
    n3 = tf.constant(1.0, dtype=tf.float64) - n1 - n2
   
    pdf_sig = n1*gauss1 + n2*gauss2 + n3*gauss3
    
    pdf_ttbar = ttbar_para_qT_cond(xvals, args[0]['parms_ttbar'])
    pdf_ewk = ewk_para_qT_cond(xvals, args[0]['parms_ewk'])
    
    qTbin = tf.cast(tf.math.multiply(xvals[0], tf.constant(2, dtype=tf.float64)), dtype=tf.int64)
    k1 = tf.gather(args[0]['fracs_ttbar'], qTbin)
    k2 = tf.gather(args[0]['fracs_ewk'], qTbin)
    k3 = tf.constant(1.0, dtype=tf.float64) - k1 - k2
    pdf = k1*pdf_ttbar + k2*pdf_ewk + k3*pdf_sig
    return pdf
    
    
    
def data_perp(xvals, p, *args):
   
    gauss1 = bf.func_gauss(xvals[0], p[3], p[0])
    gauss2 = bf.func_gauss(xvals[0], p[3], p[1])
    gauss3 = bf.func_gauss(xvals[0], p[3], p[2])
    
    n1 = tf.constant(0.1, dtype=tf.float64)
    n2 = tf.constant(0.4, dtype=tf.float64)
    n3 = tf.constant(1.0, dtype=tf.float64) - n1 - n2

    pdf_sig = n1*gauss1 + n2*gauss2 + n3*gauss3
    

    pdf_ttbar = ttbar_perp_cond([args[0]['qT'], xvals[0]], args[0]['parms_ttbar'])
    pdf_ewk = ewk_perp_cond([args[0]['qT'], xvals[0]], args[0]['parms_ewk'])

    qTbin = tf.cast(tf.math.multiply(args[0]['qT'], tf.constant(2, dtype=tf.float64)), dtype=tf.int64)
    k1 = tf.gather(args[0]['fracs_ttbar'], qTbin)
    k2 = tf.gather(args[0]['fracs_ewk'], qTbin)
    k3 = tf.constant(1.0, dtype=tf.float64) - k1 - k2
    pdf = k1*pdf_ttbar + k2*pdf_ewk + k3*pdf_sig
    
    return pdf
   
    
def data_perp_cond(xvals, p, *args):

    pdf_ttbar = ttbar_perp_cond(xvals, args[0]['parms_ttbar'])
    pdf_ewk = ewk_perp_cond(xvals, args[0]['parms_ewk'])

    sigma1 = bf.power(xvals[0], p[0], p[1], p[2])
    sigma2 = bf.power(xvals[0], p[3], p[4], 4.124451797563)
    sigma3 = bf.power(xvals[0], p[5], p[6], p[7])
 
    mean = bf.linear(xvals[0], p[8], p[9])

    
    gauss1 = bf.func_gauss(xvals[1], mean, sigma1)
    gauss2 = bf.func_gauss(xvals[1], mean, sigma2)
    gauss3 = bf.func_gauss(xvals[1], mean, sigma3)


    n1 = p[10]
    n2 = p[11]
    n3 = tf.constant(1.0, dtype=tf.float64) - n1 - n2
   
    pdf_sig = n1*gauss1 + n2*gauss2 + n3*gauss3
    

    qTbin = tf.cast(tf.math.multiply(xvals[0], tf.constant(2, dtype=tf.float64)), dtype=tf.int64)
    k1 = tf.gather(args[0]['fracs_ttbar'], qTbin)
    k2 = tf.gather(args[0]['fracs_ewk'], qTbin)
    k3 = tf.constant(1.0, dtype=tf.float64) - k1 - k2 
    pdf = k1*pdf_ttbar + k2*pdf_ewk + k3*pdf_sig
    return pdf
    

    
def data_perp_cond_mctemplate(xvals, p, *args):

    pdf_ttbar = ttbar_perp_cond(xvals, args[0]['parms_ttbar'])
    pdf_ewk = ewk_perp_cond(xvals, args[0]['parms_ewk'])
    #pdf_zz = zz_perp_cond(xvals, args[0]['parms_zz'])

    sigma1 = bf.power(xvals[0], p[0], p[1], p[2])
    sigma2 = bf.power(xvals[0], p[3], p[4], p[5])
    sigma3 = bf.power(xvals[0], p[6], p[7], p[8])
    sigma4 = bf.power(xvals[0], p[9], p[10], p[11])
    
    '''
    sigma1 = bf.power(xvals[0], 6.06060957e-01, 1.09962518e+00, 1.85775021e+01)
    sigma2 = bf.power(xvals[0], p[3], p[4], p[5])
    sigma3 = bf.power(xvals[0], p[6], p[7], p[8])
    sigma4 = bf.power(xvals[0], p[9], p[10], p[11])
    
    [ 6.06060957e-01  1.09962518e+00  1.85775021e+01  2.76467250e-01
  4.37312711e-01  3.72911531e+00  1.33359880e-01  7.08404957e-01
  7.40972012e+00  3.77204829e-02  1.05512694e+00  1.20004687e+01
 -1.14585336e-02 -1.31798159e-03  1.09951457e-03  1.10259383e-01
  6.63522155e-01]

    
    '''
    
    mean = bf.linear(xvals[0], p[12], p[13])
    #mean1 = linear(xvals[0], p[12], p[13])
    #mean2 = linear(xvals[0], p[14], p[15])
    #mean3 = linear(xvals[0], p[16], p[17])
    #mean4 = tf.constant(0, dtype=tf.float64)
    
    gauss1 = bf.func_gauss(xvals[1], mean, sigma1)
    gauss2 = bf.func_gauss(xvals[1], mean, sigma2)
    gauss3 = bf.func_gauss(xvals[1], mean, sigma3)
    gauss4 = bf.func_gauss(xvals[1], mean, sigma4)

    n1 = p[14]
    n2 = 0.21 #p[15] # 0.21 
    n3 = p[15]
    n4 = tf.constant(1.0, dtype=tf.float64) - n1 - n2 - n3
   
    pdf_sig = n1*gauss1 + n2*gauss2 + n3*gauss3 + n4*gauss4
    
    # LnN: https://www.physics.ucla.edu/~cousins/stats/cousins_lognormal_prior.pdf
    # case needed as p[x] is double during plotting?!
    #ttbar_unc = math.log(1.001)
    #ttbar_constr = func_lognormal(tf.cast(p[11], dtype=tf.float64), tf.constant(0.0, dtype=tf.float64), tf.constant(ttbar_unc, dtype=tf.float64))
    
    #ewk_unc = math.log(1.01)
    #ewk_constr = func_lognormal(tf.cast(p[12], dtype=tf.float64), tf.constant(0.0, dtype=tf.float64), tf.constant(ewk_unc, dtype=tf.float64))

    qTbin = tf.cast(tf.math.multiply(xvals[0], tf.constant(2, dtype=tf.float64)), dtype=tf.int64)
    k1 = tf.gather(args[0]['fracs_ttbar'], qTbin)
    k2 = tf.gather(args[0]['fracs_ewk'], qTbin)
    #k3 = tf.gather(args[0]['fracs_zz'], qTbin)
    k3 = tf.constant(1.0, dtype=tf.float64) - k1 - k2 
    pdf = k1*pdf_ttbar + k2*pdf_ewk + k3*pdf_sig
    return pdf
    
