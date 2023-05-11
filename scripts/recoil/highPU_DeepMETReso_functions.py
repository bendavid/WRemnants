

import math
import tensorflow as tf
import baseFunctions as bf


    

def ewk_para(xvals, p, *args):
    gauss1 = bf.func_gauss(xvals[0], p[3], p[0])
    gauss2 = bf.func_gauss(xvals[0], p[4], p[1])
    gauss3 = bf.func_gauss(xvals[0], p[5], p[2])

    n1 = tf.constant(0.03, dtype=tf.float64)
    n2 = tf.constant(0.35, dtype=tf.float64)
    n3 = tf.constant(1.0, dtype=tf.float64) - n2 - n1

    pdf = n1*gauss1 + n2*gauss2 + n3*gauss3
    return pdf
    
def ewk_para_cond(xvals, p, *args):
    sigma1 = bf.cpol4(xvals[0], p[0], p[1], p[2], p[3], p[4])
    sigma2 = bf.cpol4(xvals[0], p[5], p[6], p[7], p[8], p[9])
    sigma3 = bf.cpol4(xvals[0], p[10], p[11], p[12], p[13], p[14])
    
    mean1 = bf.cpol4(xvals[0], p[15], p[16], p[17], p[18], p[19])
    mean2 = bf.cpol4(xvals[0], p[20], p[21], p[22], p[23], p[24])
    mean3 = bf.cpol4(xvals[0], p[25], p[26], p[27], p[28], p[29])

    gauss1 = bf.func_gauss(xvals[1], mean1, sigma1)
    gauss2 = bf.func_gauss(xvals[1], mean2, sigma2)
    gauss3 = bf.func_gauss(xvals[1], mean3, sigma3)

    n1 = p[30] 
    n2 = p[31]
    n3 = tf.constant(1, dtype=tf.float64) - n2 - n1
    
    pdf = n1*gauss1 + n2*gauss2 + n3*gauss3
    return pdf
 

 
def ewk_perp(xvals, p, *args):
    gauss1 = bf.func_gauss(xvals[0], tf.constant(0, dtype=tf.float64), p[0])
    gauss2 = bf.func_gauss(xvals[0], tf.constant(0, dtype=tf.float64), p[1])
    gauss3 = bf.func_gauss(xvals[0], tf.constant(0, dtype=tf.float64), p[2])

    n1 = tf.constant(0.03, dtype=tf.float64)
    n2 = tf.constant(0.3, dtype=tf.float64)
    n3 = tf.constant(1.0, dtype=tf.float64) - n2 - n1

    pdf = n1*gauss1 + n2*gauss2 + n3*gauss3
    return pdf
 
def ewk_perp_cond(xvals, p, *args):
    sigma1 = bf.cpol2(xvals[0], p[0], p[1], p[2])
    sigma2 = bf.cpol4(xvals[0], p[3], p[4], p[5], p[6], p[7])
    sigma3 = bf.cpol4(xvals[0], p[8], p[9], p[10], p[11], p[12])

    gauss1 = bf.func_gauss(xvals[1], tf.constant(0, dtype=tf.float64), sigma1)
    gauss2 = bf.func_gauss(xvals[1], tf.constant(0, dtype=tf.float64), sigma2)
    gauss3 = bf.func_gauss(xvals[1], tf.constant(0, dtype=tf.float64), sigma3)

    n1 = p[13]
    n2 = p[14]
    n3 = tf.constant(1, dtype=tf.float64) - n2 - n1

    pdf = n1*gauss1 + n2*gauss2 + n3*gauss3
    return pdf
 
 
 
 
def ttbar_para(xvals, p, *args):
    gauss1 = bf.func_gauss(xvals[0], p[3], p[0])
    gauss2 = bf.func_gauss(xvals[0], p[4], p[1])
    gauss3 = bf.func_gauss(xvals[0], p[5], p[2])

    n1 = tf.constant(0.3, dtype=tf.float64)
    n2 = tf.constant(0.2, dtype=tf.float64)
    n3 = tf.constant(1.0, dtype=tf.float64) - n2 - n1

    pdf = n1*gauss1 + n2*gauss2 + n3*gauss3
    return pdf

def ttbar_para_cond(xvals, p, *args):
    sigma1 = bf.cpol4(xvals[0], p[0], p[1], p[2], p[3], p[4])
    sigma2 = bf.cpol4(xvals[0], p[5], p[6], p[7], p[8], p[9])
    sigma3 = bf.cpol4(xvals[0], p[10], p[11], p[12], p[13], p[14])
    
    mean1 = bf.cpol4(xvals[0], p[15], p[16], p[17], p[18], p[19])
    mean2 = bf.cpol4(xvals[0], p[20], p[21], p[22], p[23], p[24])
    mean3 = bf.cpol4(xvals[0], p[25], p[26], p[27], p[28], p[29])

    gauss1 = bf.func_gauss(xvals[1], mean1, sigma1)
    gauss2 = bf.func_gauss(xvals[1], mean2, sigma2)
    gauss3 = bf.func_gauss(xvals[1], mean3, sigma3)

    n1 = p[30] 
    n2 = p[31]
    n3 = tf.constant(1, dtype=tf.float64) - n2 - n1

    pdf = n1*gauss1 + n2*gauss2 + n3*gauss3
    return pdf

 
 
def ttbar_perp(xvals, p, *args):
    gauss1 = bf.func_gauss(xvals[0], tf.constant(0, dtype=tf.float64), p[0])
    gauss2 = bf.func_gauss(xvals[0], tf.constant(0, dtype=tf.float64), p[1])
    gauss3 = bf.func_gauss(xvals[0], tf.constant(0, dtype=tf.float64), p[2])
    n1 = tf.constant(0.3, dtype=tf.float64)
    n2 = tf.constant(0.15, dtype=tf.float64)
    n3 = tf.constant(1.0, dtype=tf.float64) - n1 - n2
    pdf = n1*gauss1 + n2*gauss2 + n3*gauss3
    return pdf
    
def ttbar_perp_cond(xvals, p, *args):
    sigma1 = bf.cpol2(xvals[0], p[0], p[1], p[2])
    sigma2 = bf.cpol2(xvals[0], p[3], p[4], p[5])
    sigma3 = bf.cpol2(xvals[0], p[6], p[7], p[8])
    gauss1 = bf.func_gauss(xvals[1], tf.constant(0, dtype=tf.float64), sigma1)
    gauss2 = bf.func_gauss(xvals[1], tf.constant(0, dtype=tf.float64), sigma2)
    gauss3 = bf.func_gauss(xvals[1], tf.constant(0, dtype=tf.float64), sigma3)
    n1 = p[9]
    n2 = p[10]
    n3 = tf.constant(1.0, dtype=tf.float64) - n1 - n2
    pdf = n1*gauss1 + n2*gauss2 + n3*gauss3
    return pdf
 

def sumet(xvals, p, *args):

    gauss1 = bf.func_gauss(xvals[0], p[3], p[0])
    gauss2 = bf.func_gauss(xvals[0], p[3], p[1])
    gauss3 = bf.func_gauss(xvals[0], p[3], p[2])
 
    
    n1 = tf.constant(0.3, dtype=tf.float64)
    n2 = tf.constant(0.3, dtype=tf.float64)
    n3 = tf.constant(1.0, dtype=tf.float64) - n1 - n2 

    pdf = n1*gauss1 + n2*gauss2 + n3*gauss3
    return pdf


def dy_para(xvals, p, *args):
    
    k0 = tf.constant(3.75252e+01, dtype=tf.float64)
    k1 = tf.constant(-5.38110e-02, dtype=tf.float64)
    k2 = tf.constant(6.33144e-01, dtype=tf.float64)
    k3 = tf.constant(-3.52892e-02, dtype=tf.float64)
    k4 = tf.constant(7.99257e-04, dtype=tf.float64)
    k5 = tf.constant(1.14259e-06, dtype=tf.float64)
    k6 = tf.constant(-2.67102e-07, dtype=tf.float64)
    k7 = tf.constant(-2.08546e-09, dtype=tf.float64)
    k8 = tf.constant(1.41637e-10, dtype=tf.float64)
    k9 = tf.constant(-9.49462e-13, dtype=tf.float64)
    mean = bf.pw_poly8_poly1(args[0]['qT'], k0, k1, k2, k3, k4, k5, k6, k7, k8, k9)
    mean = 0

    gauss1 = bf.func_gauss(xvals[0], mean+p[5], p[0])
    gauss2 = bf.func_gauss(xvals[0], mean+p[5], p[1])
    gauss3 = bf.func_gauss(xvals[0], mean+p[6], p[2])
    gauss4 = bf.func_gauss(xvals[0], mean+p[7], p[3])
    gauss5 = bf.func_gauss(xvals[0], mean+p[8], p[4])
    
    n1 = tf.constant(0.00017167274118001641, dtype=tf.float64) #  0.001
    n2 = tf.constant(0.0186479272549545, dtype=tf.float64) #  0.05
    n3 = tf.constant(0.26194731500718077, dtype=tf.float64) #  0.35
    n4 = tf.constant(0.5487736947147043, dtype=tf.float64) #  0.45
    n5 = tf.constant(1.0, dtype=tf.float64) - n1 - n2 - n3 - n4
    
    pdf = n1*gauss1 + n2*gauss2 + n3*gauss3 + n4*gauss4 + n5*gauss5
    return pdf

    
def dy_para_cond(xvals, p, *args):
    k0 = tf.constant(3.75252e+01, dtype=tf.float64)
    k1 = tf.constant(-5.38110e-02, dtype=tf.float64)
    k2 = tf.constant(6.33144e-01, dtype=tf.float64)
    k3 = tf.constant(-3.52892e-02, dtype=tf.float64)
    k4 = tf.constant(7.99257e-04, dtype=tf.float64)
    k5 = tf.constant(1.14259e-06, dtype=tf.float64)
    k6 = tf.constant(-2.67102e-07, dtype=tf.float64)
    k7 = tf.constant(-2.08546e-09, dtype=tf.float64)
    k8 = tf.constant(1.41637e-10, dtype=tf.float64)
    k9 = tf.constant(-9.49462e-13, dtype=tf.float64)
    mean = bf.pw_poly8_poly1(xvals[0], k0, k1, k2, k3, k4, k5, k6, k7, k8, k9)
    mean=0
    
    sigma1 = bf.linear(xvals[0], p[0], p[1])
    sigma2 = bf.cpol5(xvals[0], p[2], p[3], p[4], p[5], p[6], p[7])
    sigma3 = bf.cpol5(xvals[0], p[8], p[9], p[10], p[11], p[12], p[13])
    sigma4 = bf.cpol5(xvals[0], p[14], p[15], p[16], p[17], p[18], p[19])
    sigma5 = bf.cpol5(xvals[0], p[20], p[21], p[22], p[23], p[24], p[25])
    
    mean1 = bf.cpol5(xvals[0], p[26], p[27], p[28], p[29], p[30], p[31])
    mean2 = bf.cpol5(xvals[0], p[32], p[33], p[34], p[35], p[36], p[37])
    mean3 = bf.cpol5(xvals[0], p[38], p[39], p[40], p[41], p[42], p[43])
    mean4 = bf.cpol5(xvals[0], p[44], p[45], p[46], p[47], p[48], p[49])
    
    gauss1 = bf.func_gauss(xvals[1], mean1, sigma1)
    gauss2 = bf.func_gauss(xvals[1], mean1, sigma2)
    gauss3 = bf.func_gauss(xvals[1], mean2, sigma3)
    gauss4 = bf.func_gauss(xvals[1], mean3, sigma4)
    gauss5 = bf.func_gauss(xvals[1], mean4, sigma5)
    
    n1 = p[50]
    n2 = p[51]
    n3 = p[52]
    n4 = p[53]
    n5 = tf.constant(1.0, dtype=tf.float64) - n1 - n2 - n3 - n4
   
    pdf = n1*gauss1 + n2*gauss2 + n3*gauss3 + n4*gauss4 + n5*gauss5
    return pdf

 
def dy_para_cond_v1(xvals, p, *args):

    k0 = tf.constant(2.98078e+01, dtype=tf.float64)
    k1 = tf.constant(-2.50056e-02, dtype=tf.float64)
    k2 = tf.constant(5.03451e-01, dtype=tf.float64)
    k3 = tf.constant(-2.86102e-02, dtype=tf.float64)
    k4 = tf.constant(7.96316e-04, dtype=tf.float64)
    k5 = tf.constant(-8.18177e-07, dtype=tf.float64)
    k6 = tf.constant(-2.47077e-07, dtype=tf.float64)
    k7 = tf.constant(-2.73957e-09, dtype=tf.float64)
    k8 = tf.constant(1.20691e-10, dtype=tf.float64)
    mean = bf.pw_poly7_poly1(xvals[0], k0, k1, k2, k3, k4, k5, k6, k7, k8)

    sigma1 = bf.power(xvals[0], p[0], p[0]*3.6786E-02, p[0]*3.3376E+00)
    sigma2 = bf.power(xvals[0], p[1], p[1]*2.8404E+01, p[13]*5.6922E+01)
    sigma3 = bf.power(xvals[0], p[2], p[2]*1.2515E+01, p[13]*8.6612E+01)
    sigma4 = bf.power(xvals[0], p[3], p[3]*5.2845E+02, p[13]*1.1960E+02)
    
    mean1 = mean + bf.linear(xvals[0], p[4], p[5])
    mean2 = mean + bf.linear(xvals[0], p[6], p[7])
    mean3 = mean + bf.linear(xvals[0], p[8], p[9])
    mean4 = mean + bf.linear(xvals[0], p[10], p[11])
    
    gauss1 = bf.func_gauss(xvals[1], mean1, sigma1)
    gauss2 = bf.func_gauss(xvals[1], mean2, sigma2)
    gauss3 = bf.func_gauss(xvals[1], mean3, sigma3)
    gauss4 = bf.func_gauss(xvals[1], mean4, sigma4)

    n1 = p[12]
    n2 = p[13]
    n3 = p[14]
    n4 = tf.constant(1.0, dtype=tf.float64) - n1 - n2 - n3
   
    pdf = n1*gauss1 + n2*gauss2 + n3*gauss3 + n4*gauss4
    return pdf
 


def dy_para_cond_PF(xvals, p, *args):


    k0 = tf.constant( 3.62001e+01, dtype=tf.float64)
    k1 = tf.constant(-5.07213e-02, dtype=tf.float64)
    k2 = tf.constant( 6.30362e-01, dtype=tf.float64)
    k3 = tf.constant(-3.46818e-02, dtype=tf.float64)
    k4 = tf.constant( 7.56997e-04, dtype=tf.float64)
    k5 = tf.constant( 1.40997e-06, dtype=tf.float64)
    k6 = tf.constant(-2.09177e-07, dtype=tf.float64)
    k7 = tf.constant(-2.76925e-09, dtype=tf.float64)
    k8 = tf.constant( 9.26276e-11, dtype=tf.float64)
    mean = bf.pw_poly7_poly1(xvals[0], k0, k1, k2, k3, k4, k5, k6, k7, k8)

    sigma1 = bf.power(xvals[0], p[0], p[1], p[2])
    sigma2 = bf.power(xvals[0], p[3], p[4], p[5])
    sigma3 = bf.power(xvals[0], p[6], p[7], p[8])
    sigma4 = bf.power(xvals[0], p[9], p[10], p[11])
    
    mean1 = mean + bf.linear(xvals[0], p[12], p[13])
    mean2 = mean + bf.linear(xvals[0], p[14], p[15])
    mean3 = mean + bf.linear(xvals[0], p[16], p[17])
    mean4 = mean + bf.linear(xvals[0], p[18], p[19])
    
    gauss1 = bf.func_gauss(xvals[1], mean1, sigma1)
    gauss2 = bf.func_gauss(xvals[1], mean2, sigma2)
    gauss3 = bf.func_gauss(xvals[1], mean3, sigma3)
    gauss4 = bf.func_gauss(xvals[1], mean4, sigma4)

    n1 = p[20]
    n2 = p[21]
    n3 = p[22]
    n4 = tf.constant(1.0, dtype=tf.float64) - n1 - n2 - n3
   
    pdf = n1*gauss1 + n2*gauss2 + n3*gauss3 + n4*gauss4
    return pdf


 
 
def dy_perp(xvals, p, *args): 
    gauss1 = bf.func_gauss(xvals[0], 0, p[0])
    gauss2 = bf.func_gauss(xvals[0], 0, p[1])
    gauss3 = bf.func_gauss(xvals[0], 0, p[2])
    gauss4 = bf.func_gauss(xvals[0], 0, p[3])
    gauss5 = bf.func_gauss(xvals[0], 0, p[4])
    
    n1 = tf.constant(0.0005, dtype=tf.float64)
    n2 = tf.constant(0.05, dtype=tf.float64)
    n3 = tf.constant(0.40, dtype=tf.float64)
    n4 = tf.constant(0.45, dtype=tf.float64)
    n5 = tf.constant(1.0, dtype=tf.float64) - n1 - n2 - n3 - n4

    pdf = n1*gauss1 + n2*gauss2 + n3*gauss3 + n4*gauss4  + n5*gauss5
    return pdf
    
   
def dy_perp_cond(xvals, p, *args):
    sigma1 = bf.linear(xvals[0], p[0], p[1])
    sigma2 = bf.cpol5(xvals[0], p[2], p[3], p[4], p[5], p[6], p[7])
    sigma3 = bf.cpol5(xvals[0], p[8], p[9], p[10], p[11], p[12], p[13])
    sigma4 = bf.cpol5(xvals[0], p[14], p[15], p[16], p[17], p[18], p[19])
    sigma5 = bf.cpol5(xvals[0], p[20], p[21], p[22], p[23], p[24], p[25])
    
    gauss1 = bf.func_gauss(xvals[1], 0, sigma1)
    gauss2 = bf.func_gauss(xvals[1], 0, sigma2)
    gauss3 = bf.func_gauss(xvals[1], 0, sigma3)
    gauss4 = bf.func_gauss(xvals[1], 0, sigma4)
    gauss5 = bf.func_gauss(xvals[1], 0, sigma5)
    
    n1 = p[26]
    n2 = p[27]
    n3 = p[28]
    n4 = p[29]
    n5 = tf.constant(1.0, dtype=tf.float64) - n1 - n2 - n3 - n4

    pdf = n1*gauss1 + n2*gauss2 + n3*gauss3 + n4*gauss4 + n5*gauss5
    return pdf
 

def dy_perp_cond_old(xvals, p, *args):

    sigma1 = bf.power(xvals[0], p[0], p[1], p[2])
    sigma2 = bf.power(xvals[0], p[3], p[4], p[5])
    sigma3 = bf.power(xvals[0], p[6], p[7], p[8])
    sigma4 = bf.power(xvals[0], p[9], p[10], p[11])
    sigma5 = bf.power(xvals[0], p[12], p[13], p[14])
    
    gauss1 = bf.func_gauss(xvals[1], 0, sigma1)
    gauss2 = bf.func_gauss(xvals[1], 0, sigma2)
    gauss3 = bf.func_gauss(xvals[1], 0, sigma3)
    gauss4 = bf.func_gauss(xvals[1], 0, sigma4)
    gauss5 = bf.func_gauss(xvals[1], 0, sigma5)
    
    n1 = p[15]
    n2 = p[16]
    n3 = p[17]
    n4 = p[18]
    n5 = tf.constant(1.0, dtype=tf.float64) - n1 - n2 - n3 - n4
   
    pdf = n1*gauss1 + n2*gauss2 + n3*gauss3 + n4*gauss4 + n5*gauss5
    return pdf
 
 
def dy_perp_cond_v1(xvals, p, *args):

    sigma1 = bf.power(xvals[0], p[0], p[0]*(2.5287E+02), p[1])
    sigma2 = bf.power(xvals[0], p[2], p[2]*(-5.8640E+01), p[3])
    sigma3 = bf.power(xvals[0], p[4], p[4]*(3.2202E+00), p[5])
    sigma4 = bf.power(xvals[0], p[6], p[6]*(5.2319E+01), p[7])

    mean1 = bf.linear(xvals[0], p[8], p[9])
    mean2 = bf.linear(xvals[0], p[10], p[11])
    
    gauss1 = bf.func_gauss(xvals[1], mean1, sigma1)
    gauss2 = bf.func_gauss(xvals[1], mean2, sigma2)
    gauss3 = bf.func_gauss(xvals[1], mean1, sigma3)
    gauss4 = bf.func_gauss(xvals[1], mean2, sigma4)

    n1 = p[12]
    n2 = p[13]
    n3 = p[14]
    n4 = tf.constant(1.0, dtype=tf.float64) - n1 - n2 - n3
   
    pdf = n1*gauss1 + n2*gauss2 + n3*gauss3 + n4*gauss4
    return pdf
 
 

 
def data_para(xvals, p, *args):
    gauss1 = bf.func_gauss(xvals[0], p[5], p[0])
    gauss2 = bf.func_gauss(xvals[0], p[5], p[1])
    gauss3 = bf.func_gauss(xvals[0], p[6], p[2])
    gauss4 = bf.func_gauss(xvals[0], p[7], p[3])
    gauss5 = bf.func_gauss(xvals[0], p[8], p[4])
    
    n1 = tf.constant(0.0005384319600519418, dtype=tf.float64) #  0.001
    n2 = tf.constant(0.01276454399066813, dtype=tf.float64) #  0.05
    n3 = tf.constant(0.21119273740160244, dtype=tf.float64) #  0.35
    n4 = tf.constant(0.5621434052416909, dtype=tf.float64) #  0.45
    # fractions can define the crossing of sigma/mean wrt each other
    n1 = tf.constant(0.0005, dtype=tf.float64)
    n2 = tf.constant(0.05, dtype=tf.float64)
    n3 = tf.constant(0.50, dtype=tf.float64)
    n4 = tf.constant(0.35, dtype=tf.float64)
    n5 = tf.constant(1.0, dtype=tf.float64) - n1 - n2 - n3 - n4
    
    pdf_sig = n1*gauss1 + n2*gauss2 + n3*gauss3 + n4*gauss4 + n5*gauss5
    pdf_ttbar = ttbar_para_cond([args[0]['qT'], xvals[0]], args[0]['parms_ttbar'])
    pdf_ewk = ewk_para_cond([args[0]['qT'], xvals[0]], args[0]['parms_ewk'])
    
    #qTbin = tf.cast(args[0]['qT'], dtype=tf.int64) ## 
    qTbin = tf.cast(tf.math.multiply(args[0]['qT'], tf.constant(2, dtype=tf.float64)), dtype=tf.int64) # 0.5 GeV qT bins
    k1 = tf.gather(args[0]['fracs_ttbar'], qTbin)
    k2 = tf.gather(args[0]['fracs_ewk'], qTbin)
    k3 = tf.constant(1.0, dtype=tf.float64) - k1 - k2
    pdf = k1*pdf_ttbar + k2*pdf_ewk + k3*pdf_sig
    return pdf 

def data_para_cond(xvals, p, *args):
    sigma1 = bf.linear(xvals[0], p[0], p[1])
    sigma2 = bf.cpol5(xvals[0], p[2], p[3], p[4], p[5], p[6], p[7])
    sigma3 = bf.cpol5(xvals[0], p[8], p[9], p[10], p[11], p[12], p[13])
    sigma4 = bf.cpol5(xvals[0], p[14], p[15], p[16], p[17], p[18], p[19])
    sigma5 = bf.cpol5(xvals[0], p[20], p[21], p[22], p[23], p[24], p[25])
    
    mean1 = bf.cpol5(xvals[0], p[26], p[27], p[28], p[29], p[30], p[31])
    mean2 = bf.cpol5(xvals[0], p[32], p[33], p[34], p[35], p[36], p[37])
    mean3 = bf.cpol5(xvals[0], p[38], p[39], p[40], p[41], p[42], p[43])
    mean4 = bf.cpol5(xvals[0], p[44], p[45], p[46], p[47], p[48], p[49])
    
    gauss1 = bf.func_gauss(xvals[1], mean1, sigma1)
    gauss2 = bf.func_gauss(xvals[1], mean1, sigma2)
    gauss3 = bf.func_gauss(xvals[1], mean2, sigma3)
    gauss4 = bf.func_gauss(xvals[1], mean3, sigma4)
    gauss5 = bf.func_gauss(xvals[1], mean4, sigma5)
    
    n1 = p[50]
    n2 = p[51]
    n3 = p[52]
    n4 = p[53]
    n5 = tf.constant(1.0, dtype=tf.float64) - n1 - n2 - n3 - n4
    
    pdf_sig = n1*gauss1 + n2*gauss2 + n3*gauss3 + n4*gauss4 + n5*gauss5
    pdf_ttbar = ttbar_para_cond(xvals, args[0]['parms_ttbar'])
    pdf_ewk = ewk_para_cond(xvals, args[0]['parms_ewk'])
    
    #qTbin = tf.cast(xvals[0], dtype=tf.int64)
    qTbin = tf.cast(tf.math.multiply(xvals[0], tf.constant(2, dtype=tf.float64)), dtype=tf.int64) # 0.5 GeV bins
    k1 = tf.gather(args[0]['fracs_ttbar'], qTbin)
    k2 = tf.gather(args[0]['fracs_ewk'], qTbin)
    k3 = tf.constant(1.0, dtype=tf.float64) - k1 - k2
    pdf = k1*pdf_ttbar + k2*pdf_ewk + k3*pdf_sig
    return pdf
    
def data_para_cond_v1(xvals, p, *args):

    
    k0 = tf.constant(8.87682e+01, dtype=tf.float64)
    k1 = tf.constant(5.33854e-02, dtype=tf.float64)
    k2 = tf.constant(4.43228e-01, dtype=tf.float64)
    k3 = tf.constant(-1.42881e-02, dtype=tf.float64)
    k4 = tf.constant(2.77539e-04, dtype=tf.float64)
    k5 = tf.constant(-1.69169e-06, dtype=tf.float64)
    k6 = tf.constant(-2.38482e-08, dtype=tf.float64)
    k7 = tf.constant(4.20338e-10, dtype=tf.float64)
    k8 = tf.constant(-1.79513e-12, dtype=tf.float64)
    mean = bf.pw_poly7_poly1(xvals[0], k0, k1, k2, k3, k4, k5, k6, k7, k8)

    sigma1 = bf.power(xvals[0], p[0], p[0]*5.0117E-02, p[1])
    sigma2 = bf.power(xvals[0], p[2], p[2]*3.0170E+01, p[14]*7.7316E+01)
    sigma3 = bf.power(xvals[0], p[3], p[3]*2.6298E+01, p[14]*1.1472E+02)
    sigma4 = bf.power(xvals[0], p[4], p[4]*1.0389E+02, p[14]*1.5667E+02)
    
    mean1 = mean + bf.linear(xvals[0], p[5], p[6])
    mean2 = mean + bf.linear(xvals[0], p[7], p[8])
    mean3 = mean + bf.linear(xvals[0], p[9], p[10])
    mean4 = mean + bf.linear(xvals[0], p[11], p[12])
    
    gauss1 = bf.func_gauss(xvals[1], mean1, sigma1)
    gauss2 = bf.func_gauss(xvals[1], mean2, sigma2)
    gauss3 = bf.func_gauss(xvals[1], mean3, sigma3)
    gauss4 = bf.func_gauss(xvals[1], mean4, sigma4)

    n1 = p[13]
    n2 = p[14]
    n3 = p[15]
    n4 = tf.constant(1.0, dtype=tf.float64) - n1 - n2 -n3
    pdf_sig = n1*gauss1 + n2*gauss2 + n3*gauss3 + n4*gauss4
    
    pdf_ttbar = ttbar_para_cond(xvals, args[0]['parms_ttbar'])
    pdf_ewk = ewk_para_cond(xvals, args[0]['parms_ewk'])
    
    qTbin = tf.cast(xvals[0], dtype=tf.int64)
    k1 = tf.gather(args[0]['fracs_ttbar'], qTbin)
    k2 = tf.gather(args[0]['fracs_ewk'], qTbin)
    k3 = tf.constant(1.0, dtype=tf.float64) - k1 - k2
    pdf = k1*pdf_ttbar + k2*pdf_ewk + k3*pdf_sig
    return pdf
     
    
    
    
def data_perp(xvals, p, *args):
    gauss1 = bf.func_gauss(xvals[0], 0, p[0])
    gauss2 = bf.func_gauss(xvals[0], 0, p[1])
    gauss3 = bf.func_gauss(xvals[0], 0, p[2])
    gauss4 = bf.func_gauss(xvals[0], 0, p[3])
    gauss5 = bf.func_gauss(xvals[0], 0, p[4])
    
    n1 = tf.constant(0.0005, dtype=tf.float64)
    n2 = tf.constant(0.05, dtype=tf.float64)
    n3 = tf.constant(0.40, dtype=tf.float64)
    n4 = tf.constant(0.45, dtype=tf.float64)
    n5 = tf.constant(1.0, dtype=tf.float64) - n1 - n2 - n3 - n4

    pdf_sig = n1*gauss1 + n2*gauss2 + n3*gauss3 + n4*gauss4  + n5*gauss5
    pdf_ttbar = ttbar_perp_cond([args[0]['qT'], xvals[0]], args[0]['parms_ttbar'])
    pdf_ewk = ewk_perp_cond([args[0]['qT'], xvals[0]], args[0]['parms_ewk'])

    qTbin = tf.cast(args[0]['qT'], dtype=tf.int64)
    k1 = tf.gather(args[0]['fracs_ttbar'], qTbin)
    k2 = tf.gather(args[0]['fracs_ewk'], qTbin)
    k3 = tf.constant(1.0, dtype=tf.float64) - k1 - k2
    pdf = k1*pdf_ttbar + k2*pdf_ewk + k3*pdf_sig
    
    return pdf
      
def data_perp_cond(xvals, p, *args):
    sigma1 = bf.linear(xvals[0], p[0], p[1])
    sigma2 = bf.cpol5(xvals[0], p[2], p[3], p[4], p[5], p[6], p[7])
    sigma3 = bf.cpol5(xvals[0], p[8], p[9], p[10], p[11], p[12], p[13])
    sigma4 = bf.cpol5(xvals[0], p[14], p[15], p[16], p[17], p[18], p[19])
    sigma5 = bf.cpol5(xvals[0], p[20], p[21], p[22], p[23], p[24], p[25])
    
    gauss1 = bf.func_gauss(xvals[1], 0, sigma1)
    gauss2 = bf.func_gauss(xvals[1], 0, sigma2)
    gauss3 = bf.func_gauss(xvals[1], 0, sigma3)
    gauss4 = bf.func_gauss(xvals[1], 0, sigma4)
    gauss5 = bf.func_gauss(xvals[1], 0, sigma5)
    
    n1 = p[26]
    n2 = p[27]
    n3 = p[28]
    n4 = p[29]
    n5 = tf.constant(1.0, dtype=tf.float64) - n1 - n2 - n3 - n4

    pdf_sig = n1*gauss1 + n2*gauss2 + n3*gauss3 + n4*gauss4 + n5*gauss5
    pdf_ttbar = ttbar_perp_cond(xvals, args[0]['parms_ttbar'])
    pdf_ewk = ewk_perp_cond(xvals, args[0]['parms_ewk'])
    
    qTbin = tf.cast(xvals[0], dtype=tf.int64) # 1 GeV bins
    k1 = tf.gather(args[0]['fracs_ttbar'], qTbin)
    k2 = tf.gather(args[0]['fracs_ewk'], qTbin)
    k3 = tf.constant(1.0, dtype=tf.float64) - k1 - k2 
    pdf = k1*pdf_ttbar + k2*pdf_ewk + k3*pdf_sig
    return pdf
    
def data_perp_cond_v1(xvals, p, *args):

    sigma1 = bf.power(xvals[0], p[0], p[0]*1.3862E+02, p[1])
    sigma2 = bf.power(xvals[0], p[2], p[2]*(-1.8742E+01), p[10]*5.2746E+02)
    sigma3 = bf.power(xvals[0], p[3], p[3]*8.3859E+00, p[11]*2.8648E+01)
    sigma4 = bf.power(xvals[0], p[4], p[4]*3.6196E+01, p[11]*4.0329E+01)
    
    mean1 = bf.linear(xvals[0], p[5], p[6])
    mean2 = bf.linear(xvals[0], p[7], p[8])
    
    gauss1 = bf.func_gauss(xvals[1], mean1, sigma1)
    gauss2 = bf.func_gauss(xvals[1], mean2, sigma2)
    gauss3 = bf.func_gauss(xvals[1], mean1, sigma3)
    gauss4 = bf.func_gauss(xvals[1], mean2, sigma4)

    n1 = p[9]
    n2 = p[10]
    n3 = p[11]
    n4 = tf.constant(1.0, dtype=tf.float64) - n1 - n2 - n3
   
    pdf_sig = n1*gauss1 + n2*gauss2 + n3*gauss3 + n4*gauss4
    pdf_ttbar = ttbar_perp_cond(xvals, args[0]['parms_ttbar'])
    pdf_ewk = ewk_perp_cond(xvals, args[0]['parms_ewk'])
    
    qTbin = tf.cast(xvals[0], dtype=tf.int64)
    k1 = tf.gather(args[0]['fracs_ttbar'], qTbin)
    k2 = tf.gather(args[0]['fracs_ewk'], qTbin)
    k3 = tf.constant(1.0, dtype=tf.float64) - k1 - k2 
    pdf = k1*pdf_ttbar + k2*pdf_ewk + k3*pdf_sig
    return pdf
    
