

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


    
