

import math
import tensorflow as tf

def pw_poly5_power_():

    # f(x) = a + b*x + c*x^2 + d*x^3 + e*x^4 + f*x^5
    # g(x) = u*x^v + w
    fLeft, dfLeft  = "[1] + [2]*x + [3]*x*x + [4]*x*x*x + [5]*x*x*x*x + [6]*x*x*x*x*x", "[2] + 2*[3]*x + 3*[4]*x*x + 4*[5]*x*x*x + 5*[6]*x*x*x*x"
    fLeftEval, dfLeftEval = fLeft.replace("x", "[0]"), dfLeft.replace("x", "[0]")
    u = "({dfLeftEval}/([7]*TMath::Power([0], [7]-1)))".format(fLeftEval=fLeftEval, dfLeftEval=dfLeftEval)
    w = "({fLeftEval}-{u}*TMath::Power([0], [7]))".format(fLeftEval=fLeftEval, u=u)
    fRight = "{u}*TMath::Power(x, [7]) + {w}".format(u=u, w=w)
    func = "(x<[0])*1*({fLeft}) +0+ (x>[0])*1*({fRight})".format(fLeft=fLeft, fRight=fRight)
    return func   

def cte(qT, p0):
    return p0

def cte_():
    func = "[0]"
    return func
    
def power____(params):
    return params[1]*tf.math.pow(params[0], params[2]) + params[3]
    
def power(qT, p0, p1, p2):
    return p0*tf.math.pow(qT, p1) + p2
    
def powerx(qT, p0, p1, p2):
    return p0*tf.math.pow(qT, p1) + p2*qT
  
def linear(qT, p0, p1):
    return p0 + qT*p1
    
def linear0(qT, p0):
    return p0*qT
    
def linear_():
    func = "[0] + [1]*x"
    return func
    
def linear0_():
    func = "[0]*x"
    return func
    
def quadratic(qT, p0, p1, p2):
    return p0 + qT*p1 + qT*qT*p2
 
def pol5_():
    func = "[0] + [1]*x + [2]*x*x + [3]*x*x*x + [4]*x*x*x*x + [5]*x*x*x*x*x"
    return func
 
def pol6_():
    func = "[0] + [1]*x + [2]*x*x + [3]*x*x*x + [4]*x*x*x*x + [5]*x*x*x*x*x + [6]*x*x*x*x*x*x"
    return func
 
def quadratic_():
    func = "[0] + [1]*x + [2]*x*x"
    return func

def power_():
    func = "[0]*TMath::Power(x, [1]) + [2]"
    return func
    
def powerx_():
    func = "[0]*TMath::Power(x, [1]) + [2]*x"
    return func
    
def power1_():
    func = "[0]*TMath::Power(x-[2], [1])"
    return func
    
    
def pw_power_lin_():
    fLeft  = "[1]*TMath::Power(x, [2]) + [3]"
    D = "[1]*[2]*TMath::Power([0], [2]-1)"
    E = "[1]*TMath::Power([0], [2]) + [3] - [1]*[2]*TMath::Power([0], [2])"
    fRight = "((%s)*x + (%s))" % (D, E)
    func = "(x<[0])*1*({fLeft}) +0+ (x>[0])*1*({fLeft})".format(fLeft=fLeft, fRight=fRight)
    return func

    
def pw_poly4_poly1_():
    fLeft, dfLeft  = "[1] + [2]*x + [3]*x*x + [4]*x*x*x + [5]*x*x*x*x", "[2] + 2*[3]*x + 3*[4]*x*x + 4*[5]*x*x*x"
    fLeftEval, dfLeftEval = fLeft.replace("x", "[0]"), dfLeft.replace("x", "[0]")
    func = "(x<[0])*1*({fLeft}) +0+ (x>[0])*1*(  ({dfLeftEval})*x + ({fLeftEval}-({dfLeftEval})*[0]) )".format(fLeft=fLeft, fLeftEval=fLeftEval, dfLeftEval=dfLeftEval)
    return func
    
def pw_poly6_poly1_():
    fLeft, dfLeft  = "[1] + [2]*x + [3]*x*x + [4]*x*x*x + [5]*x*x*x*x + [6]*x*x*x*x*x + [7]*x*x*x*x*x*x", "[2] + 2*[3]*x + 3*[4]*x*x + 4*[5]*x*x*x + 5*[6]*x*x*x*x + 6*[7]*x*x*x*x*x"
    fLeftEval, dfLeftEval = fLeft.replace("x", "[0]"), dfLeft.replace("x", "[0]")
    func = "(x<[0])*1*({fLeft}) +0+ (x>[0])*1*(  ({dfLeftEval})*x + ({fLeftEval}-({dfLeftEval})*[0]) )".format(fLeft=fLeft, fLeftEval=fLeftEval, dfLeftEval=dfLeftEval)
    return func
    
def pol6(qT, p0, p1, p2, p3, p4, p5, p6):
    return p0 + tf.math.pow(qT, 1)*p1 + tf.math.pow(qT, 2)*p2 + tf.math.pow(qT, 3)*p3 + tf.math.pow(qT, 4)*p4 + tf.math.pow(qT, 5)*p5 + tf.math.pow(qT, 6)*p6


def heaviside(var, qT):
     
    return (tf.math.maximum(tf.constant(1.0, dtype=tf.float64), tf.math.sign(-(qT-var))))
 

def pw_poly4_poly1(qT, p0, p1, p2, p3, p4, p5):
    x = p0
    fLeft = p1 + p2*qT + p3*tf.math.pow(qT, 2) + p4*tf.math.pow(qT, 3) + p5*tf.math.pow(qT, 4)
    fLeftEval = p1 + p2*x + p3*tf.math.pow(x, 2) + p4*tf.math.pow(x, 3) + p5*tf.math.pow(x, 4)
    dfLeftEval = p2 + 2.0*p3*tf.math.pow(x, 1) + 3.0*p4*tf.math.pow(x, 2) + 4.0*p5*tf.math.pow(x, 3)
    fRight = dfLeftEval*qT + (fLeftEval-dfLeftEval)*p0
    h = tf.experimental.numpy.heaviside(qT-x, 0)
    return ((1.0-h)*fLeft + h*fRight) 
    
def pw_poly4_power(qT, p0, p1, p2, p3, p4, p5, p6):
    x = p0
    fLeft = p1 + p2*qT + p3*tf.math.pow(qT, 2) + p4*tf.math.pow(qT, 3) + p5*tf.math.pow(qT, 4)
    # fLeft = p1 + p2*qT + p3*qT*qT + p4*qT*qT*qT + p5*qT*qT*qT*qT 
    
    #return fLeft
    
    fLeftEval = p1 + p2*x + p3*tf.math.pow(x, 2) + p4*tf.math.pow(x, 3) + p5*tf.math.pow(x, 4)
    dLeftEval = p2 + 2.0*p3*x + 3.0*p4*tf.math.pow(x, 2) + 4.0*p5*tf.math.pow(x, 3)
    u = dLeftEval / (p6*tf.math.pow(x, p6-1.0))
    w = fLeftEval - u*tf.math.pow(x, p6)
    h = tf.experimental.numpy.heaviside(qT-x, 0)
    u = tf.cast(u, tf.float64)
    w = tf.cast(w, tf.float64)
    fRight = (u*tf.math.pow(qT, p6) + w)
    #return fLeft
    return ((1.0-h)*fLeft + h*fRight)
    

def pw_poly6_poly1(qT, p0, p1, p2, p3, p4, p5, p6, p7):
    x = p0
    fLeft = p1 + p2*qT + p3*tf.math.pow(qT, 2) + p4*tf.math.pow(qT, 3) + p5*tf.math.pow(qT, 4) + p6*tf.math.pow(qT, 5) + p7*tf.math.pow(qT, 6)
    fLeftEval = p1 + p2*x + p3*tf.math.pow(x, 2) + p4*tf.math.pow(x, 3) + p5*tf.math.pow(x, 4) + p6*tf.math.pow(x, 5) + p7*tf.math.pow(x, 6)
    dfLeftEval = p2 + 2.0*p3*tf.math.pow(x, 1) + 3.0*p4*tf.math.pow(x, 2) + 4.0*p5*tf.math.pow(x, 3) + 5.0*p6*tf.math.pow(x, 4) + 6.0*p7*tf.math.pow(x, 5)
    fRight = dfLeftEval*qT + (fLeftEval-dfLeftEval)*p0
    h = tf.experimental.numpy.heaviside(qT-x, 0)
    return ((1.0-h)*fLeft + h*fRight)
    
def pw_power_lin(qT, p0, p1, p2, p3):

    fLeft = p1*tf.math.pow(qT, p2) + p3
    D = p1*p2*tf.math.pow(p0, p2-1.)
    E = p1*tf.math.pow(p0, p2) + p3 - p1*p2*tf.math.pow(p0, p2)
    D = tf.cast(D, tf.float64)
    E = tf.cast(E, tf.float64)
    fRight = D*qT + E
    h = tf.experimental.numpy.heaviside(qT-p0, 0)
    return ((1.0-h)*fLeft + h*fRight)
    

   
def pw_lin(qT, p0, p1):
    fLeft = p0 + p1*qT
    return fLeft
    
def pw_quadr(qT, p0, p1, p2):
    fLeft = p0 + p1*qT + p2*tf.math.pow(qT, 2.0)
    return fLeft
    
def pw_cube(qT, p0, p1, p2, p3):
    fLeft = p0 + p1*qT + p2*tf.math.pow(qT, 2.0) + p3*tf.math.pow(qT, 3.0)
    return fLeft
    
def pw_func(qT, p0, p1, p2, p3, p4, p5, p6):
    fLeft = p0 + p1*qT + p2*tf.math.pow(qT, 2.0) + p3*tf.math.pow(qT, 3.0) + p4*tf.math.pow(qT, 4.0) + p5*tf.math.pow(qT, 7.0) + p6*tf.math.pow(qT, 6.0)
    return fLeft

def pw_poly4_power____(params):
    qT = params[0]
    x = params[7]
    fLeft = params[1] + params[2]*qT + params[3]*qT*qT + params[4]*qT*qT*qT + params[5]*qT*qT*qT*qT
    if qT < x:
        return fLeft
    else:
        fLeftEval = params[1] + params[2]*x + params[3]*x*x + params[4]*x*x*x + params[5]*x*x*x*x
        dLeftEval = params[2] + 2.0*params[3]*x + 3.0*params[4]*x*x + 4.0*params[5]*x*x*x
        u = dLeftEval / (params[6]*tf.math.pow(x, params[6]-1.0))
        w = fLeftEval - u*tf.math.pow(x, params[6])
        return (u*tf.math.pow(qT, params[6]) + w)

  
def pw_poly4_power_():
    # f(x) = a + b*x + c*x^2 + d*x^3 + e*x^4
    # g(x) = u*x^v + w
    fLeft, dfLeft  = "[1] + [2]*x + [3]*x*x + [4]*x*x*x + [5]*x*x*x*x", "[2] + 2*[3]*x + 3*[4]*x*x + 4*[5]*x*x*x"
    fLeftEval, dfLeftEval = fLeft.replace("x", "[0]"), dfLeft.replace("x", "[0]")
    A = "({dfLeftEval}/([6]*TMath::Power([0], [6]-1)))".format(fLeftEval=fLeftEval, dfLeftEval=dfLeftEval)
    #w = "({fLeftEval}-{u}*TMath::Power([0], [6]))".format(fLeftEval=fLeftEval, u=u)
    #w = "({fLeftEval}-({dfLeftEval})*[0]/[6])".format(fLeftEval=fLeftEval, dfLeftEval=dfLeftEval)
    C = "({fLeftEval}-({A})*TMath::Power([0], [6]))".format(fLeftEval=fLeftEval, A=A)
    fRight = "(  ({A})*TMath::Power(x, [6]) + {C}  )".format(A=A, C=C)
    func = "(x<[0])*1*({fLeft}) +0+ (x>[0])*1*({fRight})".format(fLeft=fLeft, fRight=fRight)
    return func
   


def func_gauss(x, mu, sigma):
    return tf.math.exp(-0.5*tf.math.pow((x-mu)/sigma, 2.))/(sigma*(2*math.pi)**0.5)

def func_lognormal(x, mu, sigma):
    return func_gauss(tf.math.log(x), mu, sigma)/x


def func_gauss1(x, mu, sigma):
    return tf.math.exp(-0.5*tf.math.pow((x-mu)/sigma, 2.))/(sigma*(2*math.pi)**0.5)


#@tf.function
def model_4gauss_cond(xvals, p):
    
    #sigma1 = pw_poly4_power(xvals[1], tf.constant(10, dtype=tf.float64), p[0], p[1], p[2], p[3], p[4], p[5])
    #sigma2 = pw_poly4_power(xvals[1], tf.constant(5, dtype=tf.float64),  p[6], p[7], p[8], p[9], p[10], p[11])
    #sigma3 = pw_poly4_power(xvals[1], tf.constant(10, dtype=tf.float64), p[12], p[13], p[14], p[15], p[16], p[17])
    #sigma4 = pw_poly4_power(xvals[1], tf.constant(5, dtype=tf.float64),  p[18], p[19], p[20], p[21], p[22], p[23])
    
    #sigma1 = pw_poly4_power(xvals[1], p[0], p[1], p[2], p[3], p[4], p[5], p[6])
    #sigma2 = pw_poly4_power(xvals[1], p[7], p[8], p[9], p[10], p[11], p[12], p[13])
    #sigma3 = pw_poly4_power(xvals[1], p[14], p[15], p[16], p[17], p[18], p[19], p[20])
    #sigma4 = pw_poly4_power(xvals[1], p[21], p[22], p[23], p[24], p[25], p[26], p[27])
    #sigma1 = pw_func(xvals[1], p[0], p[1], p[2], p[3], p[4], p[5], p[6])
    #sigma2 = pw_func(xvals[1], p[7], p[8], p[9], p[10], p[11], p[12], p[13])
    #sigma3 = pw_func(xvals[1], p[14], p[15], p[16], p[17], p[18], p[19], p[20])
    #sigma4 = pw_func(xvals[1], p[21], p[22], p[23], p[24], p[25], p[26], p[27])
    
    #sigma1 = pw_lin(xvals[1], p[0], p[1])
    #sigma2 = pw_lin(xvals[1], p[2], p[3])
    #sigma3 = pw_lin(xvals[1], p[4], p[5])
    #sigma4 = pw_lin(xvals[1], p[6], p[7])
    
    sigma1 = power(xvals[1], p[0], p[1], p[2])
    sigma2 = power(xvals[1], p[3], p[4], p[5])
    sigma3 = power(xvals[1], p[6], p[7], p[8])
    sigma4 = power(xvals[1], p[9], p[10], p[11])
    
    #sigma1 = pw_cube(xvals[1], p[0], p[1], p[2], p[3])
    #sigma2 = pw_cube(xvals[1], p[4], p[5], p[6], p[7])
    #sigma3 = pw_cube(xvals[1], p[8], p[9], p[10], p[11])
    #sigma4 = pw_cube(xvals[1], p[12], p[13], p[14], p[15])

    gauss1 = func_gauss(xvals[0], tf.constant(0, dtype=tf.float64), sigma1)
    gauss2 = func_gauss(xvals[0], tf.constant(0, dtype=tf.float64), sigma2)
    gauss3 = func_gauss(xvals[0], tf.constant(0, dtype=tf.float64), sigma3)
    gauss4 = func_gauss(xvals[0], tf.constant(0, dtype=tf.float64), sigma4)
    
    #gauss1 = func_gauss(xvals[0], 0, p[0])
    #gauss2 = func_gauss(xvals[0], 0, p[1])
    #gauss3 = func_gauss(xvals[0], 0, p[2])
    #gauss4 = func_gauss(xvals[0], 0, p[3])

    pdf = 0.08*gauss1 + 0.20*gauss2 + 0.25*gauss3 + (1.-0.08-0.2-0.25)*gauss4
    return pdf




#@tf.function
def model_4gauss_cond_norm(xvals, p):

    sigma1 = power(xvals[1], p[0], p[1], p[2])
    sigma2 = power(xvals[1], p[3], p[4], p[5])
    sigma3 = power(xvals[1], p[6], p[7], p[8])
    sigma4 = power(xvals[1], p[9], p[10], p[11])

    gauss1 = func_gauss(xvals[0], tf.constant(0, dtype=tf.float64), sigma1)
    gauss2 = func_gauss(xvals[0], tf.constant(0, dtype=tf.float64), sigma2)
    gauss3 = func_gauss(xvals[0], tf.constant(0, dtype=tf.float64), sigma3)
    gauss4 = func_gauss(xvals[0], tf.constant(0, dtype=tf.float64), sigma4)

    pdf = p[12]*gauss1 + p[13]*gauss2 + p[14]*gauss3 + (1.-p[12]-p[13]-p[14])*gauss4
    return pdf



def func_1gauss_perp(xvals, p):
    
    gauss1 = func_gauss(xvals[0], tf.constant(0, dtype=tf.float64), p[0])

    pdf = gauss1
    return pdf
    
def func_3gauss_perp(xvals, p):
    
    gauss1 = func_gauss(xvals[0], tf.constant(0, dtype=tf.float64), p[0])
    gauss2 = func_gauss(xvals[0], tf.constant(0, dtype=tf.float64), p[1])
    gauss3 = func_gauss(xvals[0], tf.constant(0, dtype=tf.float64), p[2])


    pdf = 0.7*gauss1 + 0.25*gauss2 + 0.05*gauss3
    return pdf
    
def func_2gauss_perp(xvals, p):
    
    gauss1 = func_gauss(xvals[0], tf.constant(0, dtype=tf.float64), p[0])
    gauss2 = func_gauss(xvals[0], tf.constant(0, dtype=tf.float64), p[1])

    pdf = 0.75*gauss1 + 0.25*gauss2
    return pdf
    
def func_3gauss_perp(xvals, p):
    
    gauss1 = func_gauss(xvals[0], tf.constant(0, dtype=tf.float64), p[0])
    gauss2 = func_gauss(xvals[0], tf.constant(0, dtype=tf.float64), p[1])
    gauss3 = func_gauss(xvals[0], tf.constant(0, dtype=tf.float64), p[2])

    pdf = 0.33*gauss1 + 0.33*gauss2 + (1.-0.33-0.33)*gauss3
    return pdf




def ttbar_perp(xvals, p, *args):
    
    gauss1 = func_gauss(xvals[0], tf.constant(0, dtype=tf.float64), p[0])
    gauss2 = func_gauss(xvals[0], tf.constant(0, dtype=tf.float64), p[1])
    gauss3 = func_gauss(xvals[0], tf.constant(0, dtype=tf.float64), p[2])

    n1 = tf.constant(0.5, dtype=tf.float64)
    n2 = tf.constant(0.4, dtype=tf.float64)
    n3 = tf.constant(1.0, dtype=tf.float64) - n2 - n1

    pdf = n1*gauss1 + n2*gauss2 + n3*gauss3
    return pdf
    
def ttbar_perp_cond(xvals, p, *args):

    sigma1 = quadratic(xvals[0], p[0], p[1], p[2])
    sigma2 = quadratic(xvals[0], p[3], p[4], p[5])
    sigma3 = quadratic(xvals[0], p[6], p[7], p[8])

    gauss1 = func_gauss(xvals[1], tf.constant(0, dtype=tf.float64), sigma1)
    gauss2 = func_gauss(xvals[1], tf.constant(0, dtype=tf.float64), sigma2)
    gauss3 = func_gauss(xvals[1], tf.constant(0, dtype=tf.float64), sigma3)

    n1 = p[9] #tf.constant(0.5, dtype=tf.float64)
    n2 = p[10] #tf.constant(0.49, dtype=tf.float64)
    n3 = tf.constant(1, dtype=tf.float64) - n2 - n1

    pdf = n1*gauss1 + n2*gauss2 + n3*gauss3
    return pdf
 
def ttbar_para_qT(xvals, p, *args):

    gauss1 = func_gauss(xvals[0], p[3], p[0])
    gauss2 = func_gauss(xvals[0], p[4], p[1])
    gauss3 = func_gauss(xvals[0], p[5], p[2])

    n1 = tf.constant(0.65, dtype=tf.float64)
    n2 = tf.constant(0.25, dtype=tf.float64)
    n3 = tf.constant(1.0, dtype=tf.float64) - n2 - n1

    pdf = n1*gauss1 + n2*gauss2 + n3*gauss3
    return pdf
    
    
def ttbar_para_qT_highpu(xvals, p, *args):

    gauss1 = func_gauss(xvals[0], p[3], p[0])
    gauss2 = func_gauss(xvals[0], p[4], p[1])
    gauss3 = func_gauss(xvals[0], p[3], p[2])

    n1 = tf.constant(0.65, dtype=tf.float64)
    n2 = tf.constant(0.25, dtype=tf.float64)
    n3 = tf.constant(1.0, dtype=tf.float64) - n2 - n1

    pdf = n1*gauss1 + n2*gauss2 + n3*gauss3
    return pdf
    
def ttbar_para_qT_cond_highpu(xvals, p, *args):

    sigma1 = quadratic(xvals[0], p[0], p[1], p[2])
    sigma2 = quadratic(xvals[0], p[3], p[4], p[5])
    sigma3 = quadratic(xvals[0], p[6], p[7], p[8])
    
    mean1 = quadratic(xvals[0], p[9], p[10], p[11])
    mean2 = quadratic(xvals[0], p[12], p[13], p[14])

    gauss1 = func_gauss(xvals[1], mean1, sigma1)
    gauss2 = func_gauss(xvals[1], mean2, sigma2)
    gauss3 = func_gauss(xvals[1], mean1, sigma3)

    n1 = p[15] 
    n2 = p[16]
    n3 = tf.constant(1, dtype=tf.float64) - n2 - n1

    pdf = n1*gauss1 + n2*gauss2 + n3*gauss3
    return pdf
 
      
def ttbar_para_qT_cond(xvals, p, *args):

    sigma1 = quadratic(xvals[0], p[0], p[1], p[2])
    sigma2 = quadratic(xvals[0], p[3], p[4], p[5])
    sigma3 = quadratic(xvals[0], p[6], p[7], p[8])
    
    mean1 = quadratic(xvals[0], p[9], p[10], p[11])
    mean2 = quadratic(xvals[0], p[12], p[13], p[14])
    mean3 = quadratic(xvals[0], p[15], p[16], p[17])

    gauss1 = func_gauss(xvals[1], mean1, sigma1)
    gauss2 = func_gauss(xvals[1], mean2, sigma2)
    gauss3 = func_gauss(xvals[1], mean3, sigma3)

    n1 = p[18] 
    n2 = p[19]
    n3 = tf.constant(1, dtype=tf.float64) - n2 - n1

    pdf = n1*gauss1 + n2*gauss2 + n3*gauss3
    return pdf
 
 
def ewk_para_qT(xvals, p, *args):

    mean1 = p[3]
    mean2 = p[4]
    mean3 = p[5]
    
    gauss1 = func_gauss(xvals[0], mean1, p[0])
    gauss2 = func_gauss(xvals[0], mean2, p[1])
    gauss3 = func_gauss(xvals[0], mean3, p[2])
    

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
    sigma1 = power(xvals[0], p[0], p[1], p[2])
    sigma2 = power(xvals[0], p[3], p[4], p[5])
    sigma3 = power(xvals[0], p[6], p[7], p[8])
    #sigma1 = quadratic(xvals[0], p[0], p[1], p[2])
    #sigma2 = quadratic(xvals[0], p[3], p[4], p[5])
    #sigma3 = quadratic(xvals[0], p[6], p[7], p[8])
    
    #pp = tf.constant([2.19057e+00, -2.37865e-01, 4.02425e-02, -1.34963e-03, 2.19717e-05, -1.73520e-07, 5.31624e-10], dtype=tf.float64)
    #mean = pol6(xvals[0], pp[0], pp[1], pp[2], pp[3], pp[4], pp[5], pp[6])
    #mean1 = quadratic(xvals[0], p[9], p[10], p[11])
    #mean2 = quadratic(xvals[0], p[12], p[13], p[14])
    #mean3 = quadratic(xvals[0], p[15], p[16], p[17])
    mean1 = linear(xvals[0], p[9], p[10])
    mean2 = linear(xvals[0], p[11], p[12])
    mean3 = linear(xvals[0], p[13], p[14])

    gauss1 = func_gauss(xvals[1], mean1, sigma1)
    gauss2 = func_gauss(xvals[1], mean2, sigma2)
    gauss3 = func_gauss(xvals[1], mean3, sigma3)

    #n1 = p[12]
    #n2 = p[13]
    #n3 = tf.constant(1, dtype=tf.float64) - n2 - n1
    n1 = p[15]
    n2 = p[16]
    n3 = 1.0 - n2 - n1
    
    pdf = n1*gauss1 + n2*gauss2 + n3*gauss3
    return pdf
   

def ewk_perp(xvals, p, *args):
    
    gauss1 = func_gauss(xvals[0], tf.constant(0, dtype=tf.float64), p[0])
    gauss2 = func_gauss(xvals[0], tf.constant(0, dtype=tf.float64), p[1])
    gauss3 = func_gauss(xvals[0], tf.constant(0, dtype=tf.float64), p[2])

    n1 = 0.76
    n2 = 0.22
    n3 = tf.constant(1.0, dtype=tf.float64) - n2 - n1

    pdf = n1*gauss1 + n2*gauss2 + n3*gauss3
    return pdf
 
def ewk_perp_cond(xvals, p, *args):

    sigma1 = power(xvals[0], p[0], p[1], p[2])
    sigma2 = power(xvals[0], p[3], p[4], p[5])
    sigma3 = linear(xvals[0], p[6], p[7])

    gauss1 = func_gauss(xvals[1], tf.constant(0, dtype=tf.float64), sigma1)
    gauss2 = func_gauss(xvals[1], tf.constant(0, dtype=tf.float64), sigma2)
    gauss3 = func_gauss(xvals[1], tf.constant(0, dtype=tf.float64), sigma3)

    n1 = 0.76 #p[8] #0.75 #  #
    n2 = 0.22 # p[9]# 0.20 # p[10] #
    n3 = tf.constant(1, dtype=tf.float64) - n2 - n1

    pdf = n1*gauss1 + n2*gauss2 + n3*gauss3
    return pdf
 


    
# ok
def func_4gauss_perp(xvals, p):
    
    gauss1 = func_gauss(xvals[0], tf.constant(0, dtype=tf.float64), p[0])
    gauss2 = func_gauss(xvals[0], tf.constant(0, dtype=tf.float64), p[1])
    gauss3 = func_gauss(xvals[0], tf.constant(0, dtype=tf.float64), p[2])
    gauss4 = func_gauss(xvals[0], tf.constant(0, dtype=tf.float64), p[3])

    pdf = 0.08*gauss1 + 0.20*gauss2 + 0.25*gauss3 + (1.-0.08-0.2-0.25)*gauss4
    return pdf
 




# ok   
def func_4gauss_para(xvals, p):
    
    gauss1 = func_gauss(xvals[0], p[4], p[0])
    gauss2 = func_gauss(xvals[0], p[5], p[1])
    gauss3 = func_gauss(xvals[0], p[4], p[2])
    gauss4 = func_gauss(xvals[0], p[5], p[3])
    
    n1 = 0.15 #p[8]
    n2 = 0.20 #p[9]
    n3 = 0.25 #p[10]
    n4 = 1 - n1 - n2 - n3

    #pdf = 0.1*gauss1 + 0.20*gauss2 + 0.25*gauss3 + (1.-0.1-0.2-0.25)*gauss4
    pdf = n1*gauss1 + n2*gauss2 + n3*gauss3 + n4*gauss4
    return pdf
    
    
#@tf.function
def func_4gauss_norm(xvals, p):
    
    gauss1 = func_gauss(xvals[0], 0, p[0])
    gauss2 = func_gauss(xvals[0], 0, p[1])
    gauss3 = func_gauss(xvals[0], 0, p[2])
    gauss4 = func_gauss(xvals[0], 0, p[3])

    pdf = p[4]*gauss1 + p[5]*gauss2 + p[6]*gauss3 + (1.-p[4]-p[5]-p[6])*gauss4
    return pdf
    
    
def func_3gauss_norm(xvals, p):
    
    gauss1 = func_gauss(xvals[0], 0, p[0])
    gauss2 = func_gauss(xvals[0], 0, p[1])
    gauss3 = func_gauss(xvals[0], 0, p[2])


    pdf = p[3]*gauss1 + p[4]*gauss2 + (1.-p[3]-p[4])*gauss3
    return pdf
    

# ok
def mc_para_4gauss(xvals, p):

    sigma1 = power(xvals[0], p[0], p[1], p[2])
    sigma2 = power(xvals[0], p[3], p[4], p[5])
    sigma3 = power(xvals[0], p[6], p[7], p[8])
    sigma4 = power(xvals[0], p[9], p[10], p[11])
    
    #mu1 = pw_poly4_power(xvals[0], p[12], p[13], p[14], p[15], p[16], p[17], p[18])
    #mu2 = pw_poly4_power(xvals[0], p[19], p[20], p[21], p[22], p[23], p[24], p[25])
    
    mu1 = power(xvals[0], p[12], p[13], p[14])
    mu2 = power(xvals[0], p[15], p[16], p[17])
    #mu3 = power(xvals[0], p[18], p[19], p[20])
    #mu4 = power(xvals[0], p[21], p[22], p[23])
    #mu1 = pw_power_lin(xvals[0], p[12], p[13], p[14], p[15])
    #mu2 = pw_power_lin(xvals[0], p[16], p[17], p[18], p[19])

    gauss1 = func_gauss(xvals[1], mu1, sigma1)
    gauss2 = func_gauss(xvals[1], mu2, sigma2)
    gauss3 = func_gauss(xvals[1], mu1, sigma3)
    gauss4 = func_gauss(xvals[1], mu2, sigma4)
    
    #n1 = p[18]
    #n2 = p[19]
    #n3 = p[20]
    n1 = 0.08 #p[8]
    n2 = 0.20 #p[9]
    n3 = 0.25 #p[10]
    n4 = 1 - n1 - n2 - n3

    #pdf = 0.1*gauss1 + 0.20*gauss2 + 0.25*gauss3 + (1.-0.1-0.2-0.25)*gauss4
    pdf = n1*gauss1 + n2*gauss2 + n3*gauss3 + n4*gauss4
    #pdf = p[20]*gauss1 + p[21]*gauss2 + p[22]*gauss3 + (1.-p[20]-p[21]-p[22])*gauss4
    #pdf = p[26]*gauss1 + p[27]*gauss2 + p[28]*gauss3 + (1.-p[26]-p[27]-p[28])*gauss4
    return pdf   
    
    
# ok

def mc_perp_4gauss(xvals, p):

    sigma1 = power(xvals[0], p[0], p[1], p[2])
    sigma2 = power(xvals[0], p[3], p[4], p[5])
    sigma3 = power(xvals[0], p[6], p[7], p[8])
    sigma4 = power(xvals[0], p[9], p[10], p[11])
    
    gauss1 = func_gauss(xvals[1], tf.constant(0, dtype=tf.float64), sigma1)
    gauss2 = func_gauss(xvals[1], tf.constant(0, dtype=tf.float64), sigma2)
    gauss3 = func_gauss(xvals[1], tf.constant(0, dtype=tf.float64), sigma3)
    gauss4 = func_gauss(xvals[1], tf.constant(0, dtype=tf.float64), sigma4)

    n1 = p[12]
    n2 = p[13]
    n3 = p[14]
    n4 = tf.constant(1.0, dtype=tf.float64) - n1 - n2 - n3

    pdf = n1*gauss1 + n2*gauss2 + n3*gauss3 + n4*gauss4
    return pdf
    
    
def data_perp_4gauss(xvals, p):

    sigma1 = power(xvals[0], p[0], p[1], p[2])
    sigma2 = power(xvals[0], p[3], p[4], p[5])
    sigma3 = power(xvals[0], p[6], p[7], p[8])
    sigma4 = power(xvals[0], p[9], p[10], p[11])

    gauss1 = func_gauss(xvals[1], tf.constant(0, dtype=tf.float64), sigma1)
    gauss2 = func_gauss(xvals[1], tf.constant(0, dtype=tf.float64), sigma2)
    gauss3 = func_gauss(xvals[1], tf.constant(0, dtype=tf.float64), sigma3)
    gauss4 = func_gauss(xvals[1], tf.constant(0, dtype=tf.float64), sigma4)

    n1 = p[12]
    n2 = p[13]
    n3 = p[14]
    n4 = 1. - n1 - n2 - n3

    pdf = n1*gauss1 + n2*gauss2 + n3*gauss3 + n4*gauss4
    return pdf
    
    
def data_perp_4gauss(xvals, p):

    sigma1 = power(xvals[0], p[0], p[1], p[2])
    sigma2 = power(xvals[0], p[3], p[4], p[5])
    sigma3 = power(xvals[0], p[6], p[7], p[8])
    sigma4 = power(xvals[0], p[9], p[10], p[11])

    gauss1 = func_gauss(xvals[1], tf.constant(0, dtype=tf.float64), sigma1)
    gauss2 = func_gauss(xvals[1], tf.constant(0, dtype=tf.float64), sigma2)
    gauss3 = func_gauss(xvals[1], tf.constant(0, dtype=tf.float64), sigma3)
    gauss4 = func_gauss(xvals[1], tf.constant(0, dtype=tf.float64), sigma4)

    n1 = p[12]
    n2 = p[13]
    n3 = p[14]
    n4 = 1. - n1 - n2 - n3

    pdf = n1*gauss1 + n2*gauss2 + n3*gauss3 + n4*gauss4
    return pdf
    
    
def func_4gauss_perp_data(xvals, p):
    
    gauss1 = func_gauss(xvals[0], tf.constant(0, dtype=tf.float64), p[0])
    gauss2 = func_gauss(xvals[0], tf.constant(0, dtype=tf.float64), p[1])
    gauss3 = func_gauss(xvals[0], tf.constant(0, dtype=tf.float64), p[2])
    gauss4 = func_gauss(xvals[0], tf.constant(0, dtype=tf.float64), p[3])
    
    n1 = 0.05 #p[3]
    n2 = 0.3 #p[4]
    n3 = 0.4
    n4 = 1. - n1 - n2 - n3

    pdf = n1*gauss1 + n2*gauss2 + n3*gauss3 + n4*gauss4
    return pdf
    
    
    
  

  
def dy_para_qT_highpu(xvals, p, *args):

    mean1 = p[4]
    mean2 = p[5]
    mean3 = p[6]
    mean4 = p[7]
    
    gauss1 = func_gauss(xvals[0], mean1, p[0])
    gauss2 = func_gauss(xvals[0], mean2, p[1])
    gauss3 = func_gauss(xvals[0], mean3, p[2])
    gauss4 = func_gauss(xvals[0], mean4, p[3])

    n1 = tf.constant(0.0025, dtype=tf.float64)
    n2 = tf.constant(0.03, dtype=tf.float64)
    n3 = tf.constant(0.55, dtype=tf.float64)
    n4 = tf.constant(1.0, dtype=tf.float64) - n1 - n2 - n3

    pdf = n1*gauss1 + n2*gauss2 + n3*gauss3 + n4*gauss4
    return pdf



def dy_para_qT_cond_highpu(xvals, p, *args):

    sigma1 = power(xvals[0], p[0], p[1], p[2])
    sigma2 = power(xvals[0], p[3], p[4], p[5])
    sigma3 = power(xvals[0], p[6], p[7], p[8])
    sigma4 = power(xvals[0], p[9], p[10], p[11])
    
    mean1 = power(xvals[0], p[12], p[13], p[14])
    mean2 = power(xvals[0], p[15], p[16], p[17])
    mean3 = power(xvals[0], p[18], p[19], p[20])
    mean4 = power(xvals[0], p[21], p[22], p[23])
    
    gauss1 = func_gauss(xvals[1], mean1, sigma1)
    gauss2 = func_gauss(xvals[1], mean2, sigma2)
    gauss3 = func_gauss(xvals[1], mean3, sigma3)
    gauss4 = func_gauss(xvals[1], mean4, sigma4)

    n1 = p[24]
    n2 = p[25]
    n3 = p[26]
    n4 = tf.constant(1.0, dtype=tf.float64) - n1 - n2 - n3
   
    pdf = n1*gauss1 + n2*gauss2 + n3*gauss3 + n4*gauss4
    return pdf
 
    

def dy_para_qT(xvals, p, *args):

    mean1 = p[4]
    mean2 = p[5]
    #mean3 = p[6]
    #mean4 = p[7]
    
    gauss1 = func_gauss(xvals[0], mean1, p[0])
    gauss2 = func_gauss(xvals[0], mean1, p[1])
    gauss3 = func_gauss(xvals[0], mean2, p[2])
    gauss4 = func_gauss(xvals[0], mean2, p[3])

    n1 = tf.constant(0.01, dtype=tf.float64)
    n2 = tf.constant(0.2, dtype=tf.float64)
    n3 = tf.constant(0.6, dtype=tf.float64)
    n4 = tf.constant(1.0, dtype=tf.float64) - n1 - n2 - n3

    pdf = n1*gauss1 + n2*gauss2 + n3*gauss3 + n4*gauss4
    return pdf

def dy_para_qT_DeepMETReso(xvals, p, *args):

    mean1 = p[4]
    mean2 = p[5]
    mean3 = p[6]
    mean4 = p[7]
    
    gauss1 = func_gauss(xvals[0], mean1, p[0])
    gauss2 = func_gauss(xvals[0], mean2, p[1])
    gauss3 = func_gauss(xvals[0], mean3, p[2])
    gauss4 = func_gauss(xvals[0], mean4, p[3])

    n1 = tf.constant(0.4, dtype=tf.float64)
    n2 = tf.constant(0.25, dtype=tf.float64)
    n3 = tf.constant(0.30, dtype=tf.float64)
    n4 = tf.constant(1.0, dtype=tf.float64) - n1 - n2 - n3

    pdf = n1*gauss1 + n2*gauss2 + n3*gauss3 + n4*gauss4
    return pdf

 


def dy_para_qT_cond(xvals, p, *args):

    sigma1 = power(xvals[0], p[0], p[1], p[2])
    sigma2 = power(xvals[0], p[3], p[4], p[5])
    sigma3 = power(xvals[0], p[6], p[7], p[8])
    sigma4 = power(xvals[0], p[9], p[10], p[11])
    
    mean1 = power(xvals[0], p[12], p[13], p[14])
    mean2 = power(xvals[0], p[15], p[16], p[17])
    #mean3 = power(xvals[0], p[18], p[19], p[20])
    #mean4 = power(xvals[0], p[21], p[22], p[23])
    
    gauss1 = func_gauss(xvals[1], mean1, sigma1)
    gauss2 = func_gauss(xvals[1], mean1, sigma2)
    gauss3 = func_gauss(xvals[1], mean2, sigma3)
    gauss4 = func_gauss(xvals[1], mean2, sigma4)

    n1 = p[18]
    n2 = p[19]
    n3 = p[20]
    n4 = tf.constant(1.0, dtype=tf.float64) - n1 - n2 - n3
   
    pdf = n1*gauss1 + n2*gauss2 + n3*gauss3 + n4*gauss4
    return pdf
 
 
def dy_perp_highpu(xvals, p, *args):
   
    gauss1 = func_gauss(xvals[0], p[4], p[0])
    gauss2 = func_gauss(xvals[0], p[5], p[1])
    gauss3 = func_gauss(xvals[0], p[4], p[2])
    gauss4 = func_gauss(xvals[0], p[5], p[3])

    n1 = tf.constant(0.005, dtype=tf.float64)
    n2 = tf.constant(0.2, dtype=tf.float64)
    n3 = tf.constant(0.25, dtype=tf.float64)
    n4 = tf.constant(1.0, dtype=tf.float64) - n1 - n2 - n3

    pdf = n1*gauss1 + n2*gauss2 + n3*gauss3 + n4*gauss4
    return pdf 
 
def dy_perp(xvals, p, *args):
   
    gauss1 = func_gauss(xvals[0], p[4], p[0])
    gauss2 = func_gauss(xvals[0], p[5], p[1])
    gauss3 = func_gauss(xvals[0], p[4], p[2])
    gauss4 = func_gauss(xvals[0], p[5], p[3])
    
    n1 = tf.constant(0.2, dtype=tf.float64)
    n2 = tf.constant(0.35, dtype=tf.float64)
    n3 = tf.constant(0.25, dtype=tf.float64)
    n4 = tf.constant(1.0, dtype=tf.float64) - n1 - n2 - n3

    #n1 = tf.constant(0.1, dtype=tf.float64)
    #n2 = tf.constant(0.2, dtype=tf.float64)
    #n3 = tf.constant(0.25, dtype=tf.float64)
    #n4 = tf.constant(1.0, dtype=tf.float64) - n1 - n2 - n3

    pdf = n1*gauss1 + n2*gauss2 + n3*gauss3 + n4*gauss4
    return pdf
    
    
def dy_perp_DeepMETReso(xvals, p, *args):
   
    gauss1 = func_gauss(xvals[0], p[4], p[0])
    gauss2 = func_gauss(xvals[0], p[5], p[1])
    gauss3 = func_gauss(xvals[0], p[4], p[2])
    gauss4 = func_gauss(xvals[0], p[5], p[3])


    n1 = tf.constant(0.40, dtype=tf.float64)
    n2 = tf.constant(0.25, dtype=tf.float64)
    n3 = tf.constant(0.30, dtype=tf.float64)
    n4 = tf.constant(1.0, dtype=tf.float64) - n1 - n2 - n3

    pdf = n1*gauss1 + n2*gauss2 + n3*gauss3 + n4*gauss4
    return pdf

def dy_perp_cond(xvals, p, *args):

    sigma1 = power(xvals[0], p[0], p[1], p[2])
    sigma2 = power(xvals[0], p[3], p[4], p[5])
    sigma3 = power(xvals[0], p[6], p[7], p[8])
    sigma4 = power(xvals[0], p[9], p[10], p[11])
    
    mean1 = linear(xvals[0], p[12], p[13])
    mean2 = linear(xvals[0], p[14], p[15])
    
    gauss1 = func_gauss(xvals[1], mean1, sigma1)
    gauss2 = func_gauss(xvals[1], mean2, sigma2)
    gauss3 = func_gauss(xvals[1], mean1, sigma3)
    gauss4 = func_gauss(xvals[1], mean2, sigma4)

    n1 = p[16]
    n2 = p[17]
    n3 = p[18]
    n4 = tf.constant(1.0, dtype=tf.float64) - n1 - n2 - n3
   
    pdf = n1*gauss1 + n2*gauss2 + n3*gauss3 + n4*gauss4
    return pdf
    
def dy_perp_cond_highpu(xvals, p, *args):

    sigma1 = power(xvals[0], p[0], p[1], p[2])
    sigma2 = power(xvals[0], p[3], p[4], p[5])
    sigma3 = power(xvals[0], p[6], p[7], p[8])
    sigma4 = power(xvals[0], p[9], p[10], p[11])
    
    mean1 = linear(xvals[0], p[12], p[13])
    mean2 = linear(xvals[0], p[14], p[15])
    #mean3 = linear(xvals[0], p[16], p[17])
    #mean4 = linear(xvals[0], p[18], p[19])
    
    gauss1 = func_gauss(xvals[1], mean1, sigma1)
    gauss2 = func_gauss(xvals[1], mean2, sigma2)
    gauss3 = func_gauss(xvals[1], mean1, sigma3)
    gauss4 = func_gauss(xvals[1], mean2, sigma4)

    n1 = p[16]
    n2 = p[17]
    n3 = p[18]
    n4 = tf.constant(1.0, dtype=tf.float64) - n1 - n2 - n3
   
    pdf = n1*gauss1 + n2*gauss2 + n3*gauss3 + n4*gauss4
    return pdf
    
 

def data_para_qT_highpu(xvals, p, *args):

    mean1 = p[4]
    mean2 = p[5]
    mean3 = p[6]
    mean4 = p[7]
    
    gauss1 = func_gauss(xvals[0], mean1, p[0])
    gauss2 = func_gauss(xvals[0], mean2, p[1])
    gauss3 = func_gauss(xvals[0], mean3, p[2])
    gauss4 = func_gauss(xvals[0], mean4, p[3])

    n1 = tf.constant(0.0025, dtype=tf.float64)
    n2 = tf.constant(0.03, dtype=tf.float64)
    n3 = tf.constant(0.55, dtype=tf.float64)
    n4 = tf.constant(1.0, dtype=tf.float64) - n1 - n2 - n3
    pdf_sig = n1*gauss1 + n2*gauss2 + n3*gauss3 + n4*gauss4
    
    pdf_ttbar = ttbar_para_qT_cond([args[0]['qT'], xvals[0]], args[0]['parms_ttbar'])
    pdf_ewk = ewk_para_qT_cond([args[0]['qT'], xvals[0]], args[0]['parms_ewk'])
    
    qTbin = tf.cast(tf.math.multiply(args[0]['qT'], tf.constant(2, dtype=tf.float64)), dtype=tf.int64)
    k1 = tf.gather(args[0]['fracs_ttbar'], qTbin)
    k2 = tf.gather(args[0]['fracs_ewk'], qTbin)
    k3 = tf.constant(1.0, dtype=tf.float64) - k1 - k2 
    pdf = k1*pdf_ttbar + k2*pdf_ewk  + k3*pdf_sig
    return pdf
    
def data_para_qT_cond_highpu(xvals, p, *args):

    sigma1 = power(xvals[0], p[0], p[1], p[2])
    sigma2 = power(xvals[0], p[3], p[4], p[5])
    sigma3 = power(xvals[0], p[6], p[7], p[8])
    sigma4 = power(xvals[0], p[9], p[10], p[11])
    
    mean1 = power(xvals[0], p[12], p[13], p[14])
    mean2 = power(xvals[0], p[15], p[16], p[17])
    mean3 = power(xvals[0], p[18], p[19], p[20])
    mean4 = power(xvals[0], p[21], p[22], p[23])
    
    gauss1 = func_gauss(xvals[1], mean1, sigma1)
    gauss2 = func_gauss(xvals[1], mean2, sigma2)
    gauss3 = func_gauss(xvals[1], mean3, sigma3)
    gauss4 = func_gauss(xvals[1], mean4, sigma4)

    n1 = p[24]
    n2 = p[25]
    n3 = p[26]
    n4 = tf.constant(1.0, dtype=tf.float64) - n1 - n2 - n3
    pdf_sig = n1*gauss1 + n2*gauss2 + n3*gauss3 + n4*gauss4
   
    pdf_ttbar = ttbar_para_qT_cond(xvals, args[0]['parms_ttbar'])
    pdf_ewk = ewk_para_qT_cond(xvals, args[0]['parms_ewk'])
    
    qTbin = tf.cast(tf.math.multiply(xvals[0], tf.constant(2, dtype=tf.float64)), dtype=tf.int64)
    k1 = tf.gather(args[0]['fracs_ttbar'], qTbin)
    k2 = tf.gather(args[0]['fracs_ewk'], qTbin)
    k3 = tf.constant(1.0, dtype=tf.float64) - k1 - k2
    pdf = k1*pdf_ttbar + k2*pdf_ewk + k3*pdf_sig
    return pdf
 
 
def data_para_qT(xvals, p, *args):

    pdf_ttbar = ttbar_para_cond_qT([args[0]['qT'], xvals[0]], args[0]['parms_ttbar'])
    pdf_ewk = ewk_para_qT_cond([args[0]['qT'], xvals[0]], args[0]['parms_ewk'])
    pdf_zz = zz_para_qT_cond([args[0]['qT'], xvals[0]], args[0]['parms_zz'])
    
    pp = tf.constant([2.37500e+01, 1.41561e-01, 4.30537e-01, -1.62911e-02, 4.85259e-04, -6.32757e-06, 9.75226e-01], dtype=tf.float64)
    mean = pw_poly4_power(args[0]['qT'], pp[0], pp[1], pp[2], pp[3], pp[4], pp[5], pp[6])
    mean1 = mean + p[3]
    mean2 = mean + p[4]
    mean3 = mean + p[5]
    
    
    
    gauss1 = func_gauss(xvals[0], mean1, p[0])
    gauss2 = func_gauss(xvals[0], mean2, p[1])
    gauss3 = func_gauss(xvals[0], mean3, p[2])
    #gauss4 = func_gauss(xvals[0], tf.constant(0, dtype=tf.float64), p[3])

    n1 = tf.constant(0.1, dtype=tf.float64) # 0.1
    n2 = tf.constant(0.2, dtype=tf.float64) # 0.2
    n3 = tf.constant(1.0, dtype=tf.float64) - n1 - n2
    #n3 = tf.constant(0.25, dtype=tf.float64) # 0.25
    #n4 = tf.constant(1, dtype=tf.float64) - n1 - n2 - n3
    pdf_sig = n1*gauss1 + n2*gauss2 + n3*gauss3# + n4*gauss4
    
    #k1 = args[0]['frac_ttbar']
    #k2 = args[0]['frac_ewk']
    qTbin = tf.cast(tf.math.multiply(args[0]['qT'], tf.constant(2, dtype=tf.float64)), dtype=tf.int64)
    k1 = tf.gather(args[0]['fracs_ttbar'], qTbin)
    k2 = tf.gather(args[0]['fracs_ewk'], qTbin)
    k3 = tf.gather(args[0]['fracs_zz'], qTbin)
    k4 = tf.constant(1.0, dtype=tf.float64) - k1 - k2 - k3
    pdf = k1*pdf_ttbar + k2*pdf_ewk + k3*pdf_zz + k4*pdf_sig #
    return pdf
    
def data_para_qT_cond(xvals, p, *args):

    pdf_ttbar = ttbar_para_qT_cond(xvals, args[0]['parms_ttbar'])
    pdf_ewk = ewk_para_qT_cond(xvals, args[0]['parms_ewk'])
    
    
    sigma1 = power(xvals[0], p[0], p[1], p[2])
    sigma2 = power(xvals[0], p[3], p[4], p[5])
    sigma3 = power(xvals[0], p[6], p[7], p[8])
    sigma4 = power(xvals[0], p[9], p[10], p[11])
    
    mean1 = power(xvals[0], p[12], p[13], p[14])
    mean2 = power(xvals[0], p[15], p[16], p[17])
    #mean3 = power(xvals[0], p[18], p[19], p[20])
    #mean4 = power(xvals[0], p[21], p[22], p[23])
    
    gauss1 = func_gauss(xvals[1], mean1, sigma1)
    gauss2 = func_gauss(xvals[1], mean1, sigma2)
    gauss3 = func_gauss(xvals[1], mean2, sigma3)
    gauss4 = func_gauss(xvals[1], mean2, sigma4)

    n1 = p[18]
    n2 = p[19]
    n3 = p[20]
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
    
def data_perp_highpu(xvals, p, *args):

    gauss1 = func_gauss(xvals[0], p[4], p[0])
    gauss2 = func_gauss(xvals[0], p[5], p[1])
    gauss3 = func_gauss(xvals[0], p[4], p[2])
    gauss4 = func_gauss(xvals[0], p[5], p[3])

    n1 = tf.constant(0.005, dtype=tf.float64)
    n2 = tf.constant(0.2, dtype=tf.float64)
    n3 = tf.constant(0.25, dtype=tf.float64)
    n4 = tf.constant(1.0, dtype=tf.float64) - n1 - n2 - n3

    pdf_sig = n1*gauss1 + n2*gauss2 + n3*gauss3 + n4*gauss4
    
    pdf_ttbar = ttbar_perp_cond([args[0]['qT'], xvals[0]], args[0]['parms_ttbar'])
    pdf_ewk = ewk_perp_cond([args[0]['qT'], xvals[0]], args[0]['parms_ewk'])

    qTbin = tf.cast(tf.math.multiply(args[0]['qT'], tf.constant(2, dtype=tf.float64)), dtype=tf.int64)
    k1 = tf.gather(args[0]['fracs_ttbar'], qTbin)
    k2 = tf.gather(args[0]['fracs_ewk'], qTbin)
    k3 = tf.constant(1.0, dtype=tf.float64) - k1 - k2
    pdf = k1*pdf_ttbar + k2*pdf_ewk + k3*pdf_sig
    return pdf
        
  
def data_perp(xvals, p, *args):

    gauss1 = func_gauss(xvals[0], p[4], p[0])
    gauss2 = func_gauss(xvals[0], p[5], p[1])
    gauss3 = func_gauss(xvals[0], p[6], p[2])
    gauss4 = func_gauss(xvals[0], tf.constant(0, dtype=tf.float64), p[3])

    n1 = tf.constant(0.0007841984164580838, dtype=tf.float64) # 
    n2 = tf.constant(0.21299424637789444, dtype=tf.float64)
    n3 = tf.constant(0.6160576695754821, dtype=tf.float64)
    n4 = tf.constant(1.0, dtype=tf.float64) - n1 - n2 - n3

    pdf_sig = n1*gauss1 + n2*gauss2 + n3*gauss3 + n4*gauss4
    
    pdf_ttbar = ttbar_perp_cond([args[0]['qT'], xvals[0]], args[0]['parms_ttbar'])
    pdf_ewk = ewk_perp_cond([args[0]['qT'], xvals[0]], args[0]['parms_ewk'])

    qTbin = tf.cast(tf.math.multiply(args[0]['qT'], tf.constant(2, dtype=tf.float64)), dtype=tf.int64)
    k1 = tf.gather(args[0]['fracs_ttbar'], qTbin)
    k2 = tf.gather(args[0]['fracs_ewk'], qTbin)
    k3 = tf.constant(1.0, dtype=tf.float64) - k1 - k2
    pdf = k1*pdf_ttbar + k2*pdf_ewk + k3*pdf_sig
    return pdf
    

def data_perp_cond_highpu(xvals, p, *args):

    sigma1 = power(xvals[0], p[0], p[1], p[2])
    sigma2 = power(xvals[0], p[3], p[4], p[5])
    sigma3 = power(xvals[0], p[6], p[7], p[8])
    sigma4 = power(xvals[0], p[9], p[10], p[11])
    
    mean1 = linear(xvals[0], p[12], p[13])
    mean2 = linear(xvals[0], p[14], p[15])
    
    gauss1 = func_gauss(xvals[1], mean1, sigma1)
    gauss2 = func_gauss(xvals[1], mean2, sigma2)
    gauss3 = func_gauss(xvals[1], mean1, sigma3)
    gauss4 = func_gauss(xvals[1], mean2, sigma4)

    n1 = p[16]
    n2 = p[17]
    n3 = p[18]
    n4 = tf.constant(1.0, dtype=tf.float64) - n1 - n2 - n3
   
    pdf_sig = n1*gauss1 + n2*gauss2 + n3*gauss3 + n4*gauss4
    pdf_ttbar = ttbar_perp_cond(xvals, args[0]['parms_ttbar'])
    pdf_ewk = ewk_perp_cond(xvals, args[0]['parms_ewk'])
    
    qTbin = tf.cast(tf.math.multiply(xvals[0], tf.constant(2, dtype=tf.float64)), dtype=tf.int64)
    k1 = tf.gather(args[0]['fracs_ttbar'], qTbin)
    k2 = tf.gather(args[0]['fracs_ewk'], qTbin)
    #k3 = tf.gather(args[0]['fracs_zz'], qTbin)
    k3 = tf.constant(1.0, dtype=tf.float64) - k1 - k2 
    pdf = k1*pdf_ttbar + k2*pdf_ewk + k3*pdf_sig
    return pdf
    
def data_perp_cond(xvals, p, *args):

    pdf_ttbar = ttbar_perp_cond(xvals, args[0]['parms_ttbar'])
    pdf_ewk = ewk_perp_cond(xvals, args[0]['parms_ewk'])
    #pdf_zz = zz_perp_cond(xvals, args[0]['parms_zz'])

    sigma1 = power(xvals[0], p[0], p[1], p[2])
    sigma2 = power(xvals[0], p[3], p[4], p[5])
    sigma3 = power(xvals[0], p[6], p[7], p[8])
    sigma4 = power(xvals[0], p[9], p[10], p[11])
    
    mean = linear(xvals[0], p[12], p[13])
    #mean1 = linear(xvals[0], p[12], p[13])
    #mean2 = linear(xvals[0], p[14], p[15])
    #mean3 = linear(xvals[0], p[16], p[17])
    #mean4 = tf.constant(0, dtype=tf.float64)
    
    gauss1 = func_gauss(xvals[1], mean, sigma1)
    gauss2 = func_gauss(xvals[1], mean, sigma2)
    gauss3 = func_gauss(xvals[1], mean, sigma3)
    gauss4 = func_gauss(xvals[1], mean, sigma4)

    n1 = p[14]
    n2 = 0.21 # p[15]
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
    
