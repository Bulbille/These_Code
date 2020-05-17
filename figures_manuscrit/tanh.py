#!/usr/bin/python
# -*- coding: utf8
import matplotlib.pyplot as plt 
fsize=15
plt.rc("font",size=fsize)
plt.rc("font",family='serif')

import numpy as np
from scipy.optimize import curve_fit
from numpy import linalg as LA
import scipy as sy
import scipy.special as sp
import scipy.integrate as integrate
import sys 

ls = np.linspace(0.1,1,400)
tot = np.empty(len(ls))

for i,l in enumerate(ls) :
    inte = integrate.quad(lambda x,l : (1/l*np.tanh(x/l)*(1-np.tanh(x/l)**2))**2,-np.inf,np.inf,args=(l))[0]
    tot[i] = inte

plt.plot(ls,tot)
plt.xlabel('$\\xi$')
plt.ylabel('$\sigma$')
plt.savefig('tension-superficielle.pdf')
plt.show()
