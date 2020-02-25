#!/usr/bin/python
# -*- coding: utf8

#########################################
# Modèle d'interface en histogramme
#########################################
import matplotlib.pyplot as plt
fontsize=10
plt.rc("font",size=fontsize)
plt.rc("font",family='serif')

import scipy as sy
import scipy.special as sp
import scipy.integrate as integrate
from scipy.optimize import curve_fit
import numpy as np
import numpy.linalg as linalg
from pylab import *
import os #pour les chemins de fichiers
import sys #pour l'exit et l'argument
import re

alpha   =sp.ai_zeros(3)[1]
def p(h,zeros):
    nominateur = np.power(sp.airy(abs(h)+zeros[0])[0],2)
    denominateur = 2*integrate.quad(lambda x : np.power(sp.airy(abs(x)+zeros[0])[0],2),0,np.inf)[0]
    return nominateur/denominateur

def gauss(h) :
    return np.exp(-h*h/2)/np.power(2*np.pi,0.5)

x       = np.linspace(-5,5,100)

ax = plt.subplot(111)
ax.plot(x,p(x,alpha),label="Distribution de Airy")
ax.plot(x,gauss(x),label="Distribution gaussienne")
ax.set_xlabel('z')
ax.set_ylabel('f(z)')
plt.legend()
plt.savefig('distrib.pdf')
plt.show()
