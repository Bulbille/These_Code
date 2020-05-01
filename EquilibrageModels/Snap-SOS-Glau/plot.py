#!/usr/bin/python
# -*- coding: utf8


import numpy as np
import matplotlib.pyplot as plt
fontsize= 12
plt.rc("font",size=fontsize)
plt.rc("font",family='serif')
plt.rc("text",usetex=True) #Latex
import sys
import os
import re
from scipy.optimize import curve_fit

#####################################
# Affiche l'état d'un système de 
# de taille N*M 
#####################################
dataA=np.loadtxt(sys.argv[1])
dataB=np.loadtxt(sys.argv[2])
plt.plot(dataA[:,0],dataA[:,1],label="$V(h_i)=h_i$")
plt.plot(dataB[:,0],dataB[:,1],label="$V(h_i)=-|h_i|$")

plt.xlabel('$i$')
plt.ylabel('$h_i$')
plt.legend()
plt.savefig('comp-potentiels-chimiques.pdf')
plt.show()
plt.close()

