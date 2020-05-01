#!/usr/bin/python
# -*- coding: utf8

#########################################
# Modèle d'interface en histogramme
#########################################
import matplotlib.pyplot as plt
fontsize=9
labelsize=16
plt.rc("font",size=fontsize)
plt.rc("font",family='serif')

import numpy as np
import os #pour les chemins de fichiers
import sys #pour l'exit et l'argument
import re
import scipy as sy
from scipy.optimize import curve_fit

def fitpow(x,a,b,n):
    return a*np.power(x,n)+b

data = np.loadtxt(sys.argv[1])

popt,pcov = curve_fit(fitpow,data[:,0],data[:,1],p0=[0,0,2])
print(popt)

plt.plot(data[:,0],fitpow(data[:,0],*popt),'+')

plt.plot(data[:,0],data[:,1])
plt.xlabel('$f$')
plt.ylabel('$\\frac{1}{L_X} \sum_i |h_i|^2$')
#plt.ylabel('$E(f) =\\frac{1}{L_X} \sum_i (h_i-h_{i+1})^2$')
plt.savefig('sigma-kaw-airy.pdf')
plt.show()
