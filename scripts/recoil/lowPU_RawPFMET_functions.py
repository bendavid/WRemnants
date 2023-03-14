

import math
import tensorflow as tf
import baseFunctions as bf


    

def ewk_para(xvals, p, *args):

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
    
def ewk_para_cond(xvals, p, *args):

    sigma1 = bf.power(xvals[0], p[0], p[1], p[2])
    sigma2 = bf.power(xvals[0], p[3], p[4], p[5])
    sigma3 = bf.power(xvals[0], p[6], p[7], p[8])

    mean1 = bf.linear(xvals[0], p[9], p[10])
    mean2 = bf.linear(xvals[0], p[11], p[12])
    mean3 = bf.linear(xvals[0], p[13], p[14])

    gauss1 = bf.func_gauss(xvals[1], mean1, sigma1)
    gauss2 = bf.func_gauss(xvals[1], mean2, sigma2)
    gauss3 = bf.func_gauss(xvals[1], mean3, sigma3)

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
 
 
 
 
def ttbar_para(xvals, p, *args):

    gauss1 = bf.func_gauss(xvals[0], p[3], p[0])
    gauss2 = bf.func_gauss(xvals[0], p[4], p[1])
    gauss3 = bf.func_gauss(xvals[0], p[5], p[2])

    n1 = tf.constant(0.65, dtype=tf.float64)
    n2 = tf.constant(0.25, dtype=tf.float64)
    n3 = tf.constant(1.0, dtype=tf.float64) - n2 - n1

    pdf = n1*gauss1 + n2*gauss2 + n3*gauss3
    return pdf
     
def ttbar_para_cond(xvals, p, *args):

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
 



def dy_para(xvals, p, *args):

    mean1 = p[4]
    mean2 = p[5]
    mean3 = p[6]
    
    gauss1 = bf.func_gauss(xvals[0], mean1, p[0])
    gauss2 = bf.func_gauss(xvals[0], mean2, p[1])
    gauss3 = bf.func_gauss(xvals[0], mean3, p[2])
    gauss4 = bf.func_gauss(xvals[0], mean3, p[3])
    
    n1 = tf.constant(0.1, dtype=tf.float64)
    n2 = tf.constant(0.2, dtype=tf.float64)
    n3 = tf.constant(0.25, dtype=tf.float64)
    n4 = tf.constant(1.0, dtype=tf.float64) - n1 - n2 - n3

    pdf = n1*gauss1 + n2*gauss2 + n3*gauss3 + n4*gauss4
    return pdf

def dy_para_cond(xvals, p, *args):

    sigma1 = bf.power(xvals[0], p[0], p[1], p[2])
    sigma2 = bf.power(xvals[0], p[3], p[4], p[5])
    sigma3 = bf.power(xvals[0], p[6], p[7], p[8])
    sigma4 = bf.power(xvals[0], p[9], p[10], p[11])
    
    #mean1 = 0
    mean1 = bf.linear(xvals[0], p[12], p[13])
    mean2 = bf.pw_poly1_power(xvals[0], p[14], p[15], p[16], p[17])
    mean3 = bf.pw_poly1_power(xvals[0], p[18], p[19], p[20], p[21])
    
    gauss1 = bf.func_gauss(xvals[1], mean1, sigma1)
    gauss2 = bf.func_gauss(xvals[1], mean2, sigma2)
    gauss3 = bf.func_gauss(xvals[1], mean3, sigma3)
    gauss4 = bf.func_gauss(xvals[1], mean3, sigma4)

    n1 = p[22]
    n2 = p[23]
    n3 = p[24]
    n4 = tf.constant(1.0, dtype=tf.float64) - n1 - n2 - n3
   
    pdf = n1*gauss1 + n2*gauss2 + n3*gauss3 + n4*gauss4
    return pdf
 
def dy_para_cond_v1(xvals, p, *args):

    sigma1 = bf.power(xvals[0], p[0], p[0]*(1.5579E+01), p[1])
    sigma2 = bf.power(xvals[0], p[2], p[2]*(4.8547E+00), p[3])
    sigma3 = bf.power(xvals[0], p[4], p[4]*(5.1943E+00), p[5])
    sigma4 = bf.power(xvals[0], p[6], p[6]*(1.4121E+01), p[7])
    
    mean1 = bf.linear(xvals[0], p[8], p[9])
    mean2 = bf.pw_poly1_power(xvals[0], p[10]*(6.2829E+00), p[10]*(-1.2876E-01), p[10], p[11])
    mean3 = bf.pw_poly1_power(xvals[0], p[12]*(2.7974E+00), p[12]*(-2.5097E-01), p[12], p[13])

    gauss1 = bf.func_gauss(xvals[1], mean1, sigma1)
    gauss2 = bf.func_gauss(xvals[1], mean2, sigma2)
    gauss3 = bf.func_gauss(xvals[1], mean3, sigma3)
    gauss4 = bf.func_gauss(xvals[1], mean3, sigma4)

    n1 = p[14]
    n2 = p[15]
    n3 = p[16]
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
 
def dy_perp_cond_v1(xvals, p, *args):

    sigma1 = bf.power(xvals[0], p[0], p[0]*(3.6667E-01), p[0]*(1.0716E+01))
    sigma2 = bf.power(xvals[0], p[1], p[1]*(1.2241E+00), p[2])
    sigma3 = bf.power(xvals[0], p[3], p[3]*(1.0949E+00), p[4])
    sigma4 = bf.power(xvals[0], p[5], p[5]*(1.0210E+01), p[6])
    
    mean1 = bf.linear(xvals[0], p[7], p[8])
    mean2 = bf.linear(xvals[0], p[9], p[10])
    
    gauss1 = bf.func_gauss(xvals[1], mean1, sigma1)
    gauss2 = bf.func_gauss(xvals[1], mean2, sigma2)
    gauss3 = bf.func_gauss(xvals[1], mean1, sigma3)
    gauss4 = bf.func_gauss(xvals[1], mean2, sigma4)

    n1 = p[11]
    n2 = p[12]
    n3 = p[13]
    n4 = tf.constant(1.0, dtype=tf.float64) - n1 - n2 - n3
   
    pdf = n1*gauss1 + n2*gauss2 + n3*gauss3 + n4*gauss4
    return pdf
 
 
 
def data_para(xvals, p, *args):

    gauss1 = bf.func_gauss(xvals[0], p[3], p[0])
    gauss2 = bf.func_gauss(xvals[0], p[4], p[1])
    gauss3 = bf.func_gauss(xvals[0], p[4], p[2])

    n1 = tf.constant(0.2, dtype=tf.float64)
    n2 = tf.constant(0.12, dtype=tf.float64)
    n3 = tf.constant(1.0, dtype=tf.float64) - n1 - n2
   
    pdf_sig = n1*gauss1 + n2*gauss2 + n3*gauss3
    
    pdf_ttbar = ttbar_para_cond([args[0]['qT'], xvals[0]], args[0]['parms_ttbar'])
    pdf_ewk = ewk_para_cond([args[0]['qT'], xvals[0]], args[0]['parms_ewk'])
    
    qTbin = tf.cast(tf.math.multiply(args[0]['qT'], tf.constant(2, dtype=tf.float64)), dtype=tf.int64)
    k1 = tf.gather(args[0]['fracs_ttbar'], qTbin)
    k2 = tf.gather(args[0]['fracs_ewk'], qTbin)
    k3 = tf.constant(1.0, dtype=tf.float64) - k1 - k2
    pdf = k1*pdf_ttbar + k2*pdf_ewk + k3*pdf_sig
    return pdf
 
def data_para_cond(xvals, p, *args):

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

    pdf_ttbar = ttbar_para_cond(xvals, args[0]['parms_ttbar'])
    pdf_ewk = ewk_para_cond(xvals, args[0]['parms_ewk'])
    
    qTbin = tf.cast(tf.math.multiply(xvals[0], tf.constant(2, dtype=tf.float64)), dtype=tf.int64)
    k1 = tf.gather(args[0]['fracs_ttbar'], qTbin)
    k2 = tf.gather(args[0]['fracs_ewk'], qTbin)
    k3 = tf.constant(1.0, dtype=tf.float64) - k1 - k2
    pdf = k1*pdf_ttbar + k2*pdf_ewk + k3*pdf_sig
    return pdf
    
def data_para_cond_v1(xvals, p, *args):
 
    
    sigma1 = bf.power(xvals[0], p[0], p[0]*(2.2261E+02), p[1])
    sigma2 = bf.power(xvals[0], p[2], p[2]*(9.4597E+00), p[3])
    sigma3 = bf.power(xvals[0], p[4], p[4]*(3.0800E+00), p[5])
    
    mean1 = bf.power(xvals[0], p[6], p[6]*(6.8832E+00), p[6]*(1.0499E+01))
    mean2 = bf.power(xvals[0], p[7], p[7]*(2.6674E-01), p[7]*(-8.4509E-01))

    gauss1 = bf.func_gauss(xvals[1], mean1, sigma1)
    gauss2 = bf.func_gauss(xvals[1], mean2, sigma2)
    gauss3 = bf.func_gauss(xvals[1], mean2, sigma3)

    n1 = p[8]
    n2 = p[9]
    n3 = tf.constant(1.0, dtype=tf.float64) - n1 - n2
   
    pdf_sig = n1*gauss1 + n2*gauss2 + n3*gauss3
    
    pdf_ttbar = ttbar_para_cond(xvals, args[0]['parms_ttbar'])
    pdf_ewk = ewk_para_cond(xvals, args[0]['parms_ewk'])
    
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
    
    n1 = tf.constant(0.15, dtype=tf.float64)
    n2 = tf.constant(0.45, dtype=tf.float64)
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
    sigma2 = bf.power(xvals[0], p[3], p[4], p[5])
    sigma3 = bf.power(xvals[0], p[6], p[7], p[8])
    mean = bf.linear(xvals[0], p[9], p[10])

    gauss1 = bf.func_gauss(xvals[1], mean, sigma1)
    gauss2 = bf.func_gauss(xvals[1], mean, sigma2)
    gauss3 = bf.func_gauss(xvals[1], mean, sigma3)

    n1 = p[11]
    n2 = p[12]
    n3 = tf.constant(1.0, dtype=tf.float64) - n1 - n2
    
    pdf_sig = n1*gauss1 + n2*gauss2 + n3*gauss3

    


    qTbin = tf.cast(tf.math.multiply(xvals[0], tf.constant(2, dtype=tf.float64)), dtype=tf.int64)
    k1 = tf.gather(args[0]['fracs_ttbar'], qTbin)
    k2 = tf.gather(args[0]['fracs_ewk'], qTbin)
    k3 = tf.constant(1.0, dtype=tf.float64) - k1 - k2 
    pdf = k1*pdf_ttbar + k2*pdf_ewk + k3*pdf_sig
    return pdf
    
def data_perp_cond_v1(xvals, p, *args):

    sigma1 = bf.power(xvals[0], p[0], p[0]*2.8216E+01, p[1])
    sigma2 = bf.power(xvals[0], p[2], p[2]*7.4621E+00, p[3])
    sigma3 = bf.power(xvals[0], p[4], p[4]*5.2991E+00, p[5])
    mean = bf.linear(xvals[0], p[6], p[7])

    gauss1 = bf.func_gauss(xvals[1], mean, sigma1)
    gauss2 = bf.func_gauss(xvals[1], mean, sigma2)
    gauss3 = bf.func_gauss(xvals[1], mean, sigma3)


    n1 = p[8]
    n2 = p[9]
    n3 = tf.constant(1.0, dtype=tf.float64) - n1 - n2
    
    pdf_sig = n1*gauss1 + n2*gauss2 + n3*gauss3

    pdf_ttbar = ttbar_perp_cond(xvals, args[0]['parms_ttbar'])
    pdf_ewk = ewk_perp_cond(xvals, args[0]['parms_ewk'])


    qTbin = tf.cast(tf.math.multiply(xvals[0], tf.constant(2, dtype=tf.float64)), dtype=tf.int64)
    k1 = tf.gather(args[0]['fracs_ttbar'], qTbin)
    k2 = tf.gather(args[0]['fracs_ewk'], qTbin)
    k3 = tf.constant(1.0, dtype=tf.float64) - k1 - k2 
    pdf = k1*pdf_ttbar + k2*pdf_ewk + k3*pdf_sig
    return pdf
    

    
    
