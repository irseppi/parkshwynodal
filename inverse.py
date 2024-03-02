import numpy as np
from numpy.linalg import inv
#from sympy import symbols, diff
def invert_f(f0,v0,l,tprime0,fobs):
    c = 343 # m/sec speed of sound
    tprime = np.sqrt(l**2-(v0*t)**2/c)
    t = (tprime - np.sqrt(tprime**2-(1-v0**2/c**2)(tprime**2-l**2/c**2)))/(1-v0**2/c**2)

    f = f0*1/(1+(v0/c)*(v0*t/(np.sqrt(l**2-(v0*t)**2))))

    m0 = [f0, v0, l, tprime0].T


    #derivative with respect to f0
    f_prime = np.diff(f, f0)
    f_derivef0 = 1/(1+(v0/c)*(v0*t/(np.sqrt(l**2-(v0*t)**2))))
    #derivative of f with respect to v0
    f_prime = np.diff(f, v0)
    f_derivev0 = c*f0*t*v0*(t**2*v0**2-l**2)*1/(np.sqrt(l**2-(v0*t)**2)*(c*np.sqrt(l**2-(v0*t)**2)+t*v0**2)**2)
    #derivative of f with respect to l
    f_prime = np.diff(f, l)
    f_derivel = (c*f0*t*v0**2*l)*1/(np.sqrt(l**2-(v0*t)**2)*(c*np.sqrt(l**2-(v0*t)**2)+t*v0**2)**2)
    #derivative of f with respect to tprime
    f_prime = np.diff(f, tprime)
    f_derivetprime = (c**3*f0*l**2*v0**2*(c*np.sqrt((-(v0**2*tprime**2)*c**2)+tprime**2-(l**2/c**2)))+(v0**2-c**2)*tprime)/(np.sqrt((-(v0**2*tprime**2)*c**2)+tprime**2-(l**2/c**2))*np.sqrt(l**2-v0*(-np.sqrt(((-(v0**2*tprime**2)/c**2)+tprime**2-(l**2/c**2)))+tprime-(v0**2/c**2))**2)*(c**3*np.sqrt(l**2-v0*(-np.sqrt(((-(v0**2*tprime**2)/c**2)+tprime**2-(l**2/c**2)))+tprime-(v0**2/c**2)))-c**2*v0**2*np.sqrt((-(v0**2*tprime**2)/c**2)+tprime**2-(l**2/c**2))+c**2*v0**2*tprime-v0**4)**2)

    #partial derivative matrix of f with respect to m when m=m0
    G = np.array([[f_derivef0, 0, 0, 0], 
                [0, f_derivev0, 0, 0], 
                [0, 0, f_derivel, 0], 
                [0, 0, 0, f_derivetprime]])
    m = m0 +inv(G.T@G)@G.T[fobs(tprime) - f(m0,tprime)]