#!/usr/bin/python
# -*- coding: utf8


import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
import sys
import os
import re

#####################################
# Affiche l'état d'un système de 
# de taille N*M 
#####################################
def expo(x,a,b):
    return a*np.power(x,b)

data = np.loadtxt(sys.argv[1])

popt,pcov = curve_fit(expo,data[:,0],data[:,1])


coef = '{:0.2f}'.format(popt[0])
exposant = '{:0.2f}'.format(popt[1])

plt.plot(data[:,0],data[:,1],label="L = 4096")
plt.plot(data[:,0],expo(data[:,0],*popt),label='$'+coef+'t^{'+exposant+'}$')
plt.xscale('log')
plt.yscale('log')
plt.xlabel('Time in Monte Carlo steps')
plt.ylabel('Mean width $w = \sqrt{<h^2>-<h>^2}$')
plt.legend()

plt.savefig('glau-width-time.pdf')
plt.show()
plt.close()
