

import math
import tensorflow as tf

# Chebyshev polynomials 
def cpol2(qT, p0, p1, p2):
    min_ = 0
    max_ = 100
    qT_ = ((qT-min_)-(max_-qT))/(max_-min_)
    return p0 + p1*qT_ + p2*(2.*tf.math.pow(qT_,2)-1.)
    
def cpol2_():
    func = "[0] + [1]*x + [2]*(2.*x*x-1.)"
    return func

def cpol3(qT, p0, p1, p2, p3):
    min_ = 0
    max_ = 100
    qT_ = ((qT-min_)-(max_-qT))/(max_-min_)
    return p0 + p1*qT_ + p2*(2.*tf.math.pow(qT_,2)-1.) + p3*(4.*tf.math.pow(qT_,3)-3.*qT_)
    
def cpol3_():
    func = "[0] + [1]*x + [2]*(2.*x*x-1.) + [3]*(4.*x*x*x-3*x)"
    return func

def cpol4(qT, p0, p1, p2, p3, p4):
    min_ = 0
    max_ = 100
    qT_ = ((qT-min_)-(max_-qT))/(max_-min_)
    return p0 + p1*qT_ + p2*(2.*tf.math.pow(qT_,2)-1.) + p3*(4.*tf.math.pow(qT_,3)-3.*qT_) + p4*(8.*tf.math.pow(qT_,4)-8.*tf.math.pow(qT_,2)+1.)
    
def cpol4_():
    func = "[0] + [1]*x + [2]*(2.*x*x-1.) + [3]*(4.*x*x*x-3*x) + [4]*(8.*x*x*x*x-8.*x*x+1.)"
    return func

def cpol5(qT, p0, p1, p2, p3, p4, p5):
    min_ = 0
    max_ = 100
    qT_ = ((qT-min_)-(max_-qT))/(max_-min_)
    return p0 + p1*qT_ + p2*(2.*tf.math.pow(qT_,2)-1.) + p3*(4.*tf.math.pow(qT_,3)-3.*qT_) + p4*(8.*tf.math.pow(qT_,4)-8.*tf.math.pow(qT_,2)+1.) + p5*(16.*tf.math.pow(qT_,5)-20.*tf.math.pow(qT_,3)+5.*qT_)
    
def cpol5_():
    func = "[0] + [1]*x + [2]*(2.*x*x-1.) + [3]*(4.*x*x*x-3*x) + [4]*(8.*x*x*x*x-8.*x*x+1.) + [5]*(16.*x*x*x*x*x-20.*x*x*x+5.*x)"
    return func


def cpol6(qT, p0, p1, p2, p3, p4, p5, p6):
    min_ = 0
    max_ = 100
    qT_ = ((qT-min_)-(max_-qT))/(max_-min_)
    return p0 + p1*qT_ + p2*(2.*tf.math.pow(qT_,2)-1.) + p3*(4.*tf.math.pow(qT_,3)-3.*qT_) + p4*(8.*tf.math.pow(qT_,4)-8.*tf.math.pow(qT_,2)+1.) + p5*(16.*tf.math.pow(qT_,5)-20.*tf.math.pow(qT_,3)+5.*qT_) + p6*(32.*tf.math.pow(qT_,6)-48.*tf.math.pow(qT_,4)+18.*tf.math.pow(qT_,2)-1.)
    
def cpol6_():
    func = "[0] + [1]*x + [2]*(2.*x*x-1.) + [3]*(4.*x*x*x-3*x) + [4]*(8.*x*x*x*x-8.*x*x+1.) + [5]*(16.*x*x*x*x*x-20.*x*x*x+5.*x) + [6]*(32.*x*x*x*x*x*x-48.*x*x*x*x+18.*x*x-1.)"
    return func


def pw_poly1_power_():
    # f(x) = a + b*x + c*x^2 + d*x^3 + e*x^4
    # g(x) = u*x^v + w
    #fLeft, dfLeft  = "[1] + [2]*x + [3]*x*x + [4]*x*x*x + [5]*x*x*x*x", "[2] + 2*[3]*x + 3*[4]*x*x + 4*[5]*x*x*x"
    fLeft, dfLeft  = "[1] + [2]*x", "[2]"
    fLeftEval, dfLeftEval = fLeft.replace("x", "[0]"), dfLeft.replace("x", "[0]")
    A = "({dfLeftEval}/([3]*TMath::Power([0], [3]-1)))".format(fLeftEval=fLeftEval, dfLeftEval=dfLeftEval)
    #w = "({fLeftEval}-{u}*TMath::Power([0], [6]))".format(fLeftEval=fLeftEval, u=u)
    #w = "({fLeftEval}-({dfLeftEval})*[0]/[6])".format(fLeftEval=fLeftEval, dfLeftEval=dfLeftEval)
    C = "({fLeftEval}-({A})*TMath::Power([0], [3]))".format(fLeftEval=fLeftEval, A=A)
    fRight = "(  ({A})*TMath::Power(x, [3]) + {C}  )".format(A=A, C=C)
    func = "(x<[0])*1*({fLeft}) +0+ (x>[0])*1*({fRight})".format(fLeft=fLeft, fRight=fRight)
    return func


def pw_poly1_power(qT, p0, p1, p2, p3):
    x = p0
    fLeft = p1 + p2*qT
    
    fLeftEval = p1 + p2*x
    dLeftEval = p2
    u = dLeftEval / (p3*tf.math.pow(x, p3-1.0))
    w = fLeftEval - u*tf.math.pow(x, p3)
    h = tf.experimental.numpy.heaviside(qT-x, 0)
    u = tf.cast(u, tf.float64)
    w = tf.cast(w, tf.float64)
    fRight = (u*tf.math.pow(qT, p3) + w)
    #return fLeft
    return ((1.0-h)*fLeft + h*fRight)
    
    

def pw_poly4_poly1(qT, p0, p1, p2, p3, p4, p5):
    x = p0
    fLeft = p1 + p2*qT + p3*tf.math.pow(qT, 2) + p4*tf.math.pow(qT, 3) + p5*tf.math.pow(qT, 4)
    fLeftEval = p1 + p2*x + p3*tf.math.pow(x, 2) + p4*tf.math.pow(x, 3) + p5*tf.math.pow(x, 4)
    dfLeftEval = p2 + 2.0*p3*tf.math.pow(x, 1) + 3.0*p4*tf.math.pow(x, 2) + 4.0*p5*tf.math.pow(x, 3)
    fRight = dfLeftEval*qT + (fLeftEval-dfLeftEval)*p0
    h = tf.experimental.numpy.heaviside(qT-x, 0)
    return ((1.0-h)*fLeft + h*fRight) 

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
    
def pol4(qT, p0, p1, p2, p3, p4):
    return p0 + qT*p1 + qT*qT*p2 + qT*qT*qT*p3 + qT*qT*qT*qT*p4
    

    
def pol4_():
    func = "[0] + [1]*x + [2]*x*x + [3]*x*x*x + [4]*x*x*x*x"
    return func
    

 
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
    
def powerr_():
    func = "[0]*TMath::Power(x-[1], [2]) + [3]"
    return func
    
def powerx_():
    func = "[0]*TMath::Power(x, [1]) + [2]*x"
    return func
    
def power1_():
    func = "[0]*TMath::Power(x-[2], [1])"
    return func
    
def log_():
    func = "TMath::Log([0]*x-[1]) + [2]"
    return func

def log(qT, p0, p1, p2):
    return p0*tf.math.log(qT-p1) + p2
    
def pw_quadr_lin_():
    val="20"
    fLeft, dfLeft = "[0] + [1]*x + [2]*x*x", "[1] + 2*[2]*x"
    fLeftEval, dfLeftEval = fLeft.replace("x", val), dfLeft.replace("x", val)
    fRight = f"({dfLeftEval})*(x-{val}) + ({fLeftEval})"
    
    func = f"(x<{val})*1*({fLeft}) +0+ (x>{val})*1*({fRight})".format(fLeft=fLeft, fRight=fRight)
    return func
    

def pw_quadr_lin(qT, p0, p1, p2, p3):
    fLeft = p0 + p1*qT + p2*qT*qT
    fLeftEval = p0 + p1*p3 + p2*p3*p3
    dfLeftEval = p1 + tf.constant(2.0, dtype=tf.float64)*p2*p3
    dfLeftEval = tf.cast(dfLeftEval, tf.float64)
    fRight = dfLeftEval*(qT-p3) + fLeftEval
    
    fLeft = tf.cast(fLeft, tf.float64)
    fRight = tf.cast(fRight, tf.float64)
    h = tf.experimental.numpy.heaviside(qT-p3, 0)
    return ((tf.constant(1.0, dtype=tf.float64)-h)*fLeft + h*fRight)   
    

    
def potential_():

    func = "TMath::Power(x+[0], [1]) + 1/(x+[0]) + [3]"
    return func

def rational(qT, p0, p1, p2, p3):
    
    ret = p0*p1*tf.math.pow(qT+p1, p2)/(qT+p3)
    return ret
    
    
    
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
    
def pw_poly7_poly1_():
    fLeft, dfLeft  = "[1] + [2]*x + [3]*x*x + [4]*x*x*x + [5]*x*x*x*x + [6]*x*x*x*x*x + [7]*x*x*x*x*x*x + [8]*x*x*x*x*x*x*x", "[2] + 2*[3]*x + 3*[4]*x*x + 4*[5]*x*x*x + 5*[6]*x*x*x*x + 6*[7]*x*x*x*x*x + 7*[8]*x*x*x*x*x*x"
    fLeftEval, dfLeftEval = fLeft.replace("x", "[0]"), dfLeft.replace("x", "[0]")
    func = "(x<[0])*1*({fLeft}) +0+ (x>[0])*1*(  ({dfLeftEval})*x + ({fLeftEval}-({dfLeftEval})*[0]) )".format(fLeft=fLeft, fLeftEval=fLeftEval, dfLeftEval=dfLeftEval)
    return func
    
def pw_poly8_poly1_():
    fLeft, dfLeft  = "[1] + [2]*x + [3]*x*x + [4]*x*x*x + [5]*x*x*x*x + [6]*x*x*x*x*x + [7]*x*x*x*x*x*x + [8]*x*x*x*x*x*x*x + [9]*x*x*x*x*x*x*x*x", "[2] + 2*[3]*x + 3*[4]*x*x + 4*[5]*x*x*x + 5*[6]*x*x*x*x + 6*[7]*x*x*x*x*x + 7*[8]*x*x*x*x*x*x + 8*[9]*x*x*x*x*x*x*x"
    fLeftEval, dfLeftEval = fLeft.replace("x", "[0]"), dfLeft.replace("x", "[0]")
    func = "(x<[0])*1*({fLeft}) +0+ (x>[0])*1*(  ({dfLeftEval})*x + ({fLeftEval}-({dfLeftEval})*[0]) )".format(fLeft=fLeft, fLeftEval=fLeftEval, dfLeftEval=dfLeftEval)
    return func
    
def quadr_power_():
    fLeft, dfLeft  = "[1] + [2]*x + [3]*x*x", "[2] + 2*[3]*x"
    fLeftEval, dfLeftEval = fLeft.replace("x", "[0]"), dfLeft.replace("x", "[0]")
    A = f"({dfLeftEval})/([4]*TMath::Power([0], [4]-1))"
    C = f"({fLeft})-({A})*TMath::Power([0], [4])"
    g = f"({A})*TMath::Power(x, [4]) + {C}"
    func = "(x<[0])*1*({fLeft}) +0+ (x>[0])*1*({g})".format(fLeft=fLeft, fLeftEval=fLeftEval, dfLeftEval=dfLeftEval, g=g)
    return func
    
def quadr_power_cte_():
    val="10"
    fLeft, dfLeft  = "[0] + [1]*x + [2]*x*x", "[1] + 2*[2]*x"
    fLeftEval, dfLeftEval = fLeft.replace("x", val), dfLeft.replace("x", val)
    A = f"({dfLeftEval})/([3]*TMath::Power({val}, [3]-1))"
    C = f"({fLeft})-({A})*TMath::Power({val}, [3])"
    g = f"({A})*TMath::Power(x, [3]) + {C}"
    func = "(x<{val})*1*({fLeft}) +0+ (x>{val})*1*({g})".format(fLeft=fLeft, fLeftEval=fLeftEval, dfLeftEval=dfLeftEval, g=g, val=val)
    return func

def quadr_power(qT, p0, p1, p2, p3, p4):
    x = p0
    
    fLeft = p1 + p2*qT + p3*tf.math.pow(qT, 2)
    fLeftEval = p1 + p2*x + p3*tf.math.pow(x, 2)
    dLeftEval = p2 + 2.0*p3*x
    
    A = dLeftEval / (p4*tf.math.pow(x, p4-1.0))
    C = fLeftEval - A*tf.math.pow(x, p4)
    h = tf.experimental.numpy.heaviside(qT-x, 0)
    A = tf.cast(A, tf.float64)
    C = tf.cast(C, tf.float64)
    fRight = (A*tf.math.pow(qT, p4) + C)
    return ((1.0-h)*fLeft + h*fRight)
    
def quadr_power_cte(qT, p0, p1, p2, p3):
    x = tf.constant(10.0, dtype=tf.float64)
    
    fLeft = p0 + p1*qT + p2*tf.math.pow(qT, 2)
    fLeftEval = p0 + p1*x + p2*tf.math.pow(x, 2)
    dLeftEval = p1 + 2.0*p2*x
    
    A = dLeftEval / (p3*tf.math.pow(x, p3-1.0))
    C = fLeftEval - A*tf.math.pow(x, p3)
    h = tf.experimental.numpy.heaviside(qT-x, 0)
    A = tf.cast(A, tf.float64)
    C = tf.cast(C, tf.float64)
    fRight = (A*tf.math.pow(qT, p3) + C)
    return ((1.0-h)*fLeft + h*fRight)
    
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
    

def pw_poly7_poly1(qT, p0, p1, p2, p3, p4, p5, p6, p7, p8):
    x = p0
    fLeft = p1 + p2*qT + p3*tf.math.pow(qT, 2) + p4*tf.math.pow(qT, 3) + p5*tf.math.pow(qT, 4) + p6*tf.math.pow(qT, 5) + p7*tf.math.pow(qT, 6) + p8*tf.math.pow(qT, 7)
    fLeftEval = p1 + p2*x + p3*tf.math.pow(x, 2) + p4*tf.math.pow(x, 3) + p5*tf.math.pow(x, 4) + p6*tf.math.pow(x, 5) + p7*tf.math.pow(x, 6) + p8*tf.math.pow(x, 7)
    dfLeftEval = p2 + 2.0*p3*tf.math.pow(x, 1) + 3.0*p4*tf.math.pow(x, 2) + 4.0*p5*tf.math.pow(x, 3) + 5.0*p6*tf.math.pow(x, 4) + 6.0*p7*tf.math.pow(x, 5) + 7.0*p8*tf.math.pow(x, 6)
    fRight = dfLeftEval*qT + fLeftEval-dfLeftEval*p0
    h = tf.experimental.numpy.heaviside(qT-x, 0)
    return ((1.0-h)*fLeft + h*fRight)


def pw_poly8_poly1(qT, p0, p1, p2, p3, p4, p5, p6, p7, p8, p9):
    x = p0
    fLeft = p1 + p2*qT + p3*tf.math.pow(qT, 2) + p4*tf.math.pow(qT, 3) + p5*tf.math.pow(qT, 4) + p6*tf.math.pow(qT, 5) + p7*tf.math.pow(qT, 6) + p8*tf.math.pow(qT, 7) + p9*tf.math.pow(qT, 8)
    fLeftEval = p1 + p2*x + p3*tf.math.pow(x, 2) + p4*tf.math.pow(x, 3) + p5*tf.math.pow(x, 4) + p6*tf.math.pow(x, 5) + p7*tf.math.pow(x, 6) + p8*tf.math.pow(x, 7) + p9*tf.math.pow(x, 8)
    dfLeftEval = p2 + 2.0*p3*tf.math.pow(x, 1) + 3.0*p4*tf.math.pow(x, 2) + 4.0*p5*tf.math.pow(x, 3) + 5.0*p6*tf.math.pow(x, 4) + 6.0*p7*tf.math.pow(x, 5) + 7.0*p8*tf.math.pow(x, 6) + 8.0*p9*tf.math.pow(x, 7)
    fRight = dfLeftEval*qT + fLeftEval-dfLeftEval*p0
    h = tf.experimental.numpy.heaviside(qT-x, 0)
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


    
