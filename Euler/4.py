#!/usr/bin/python
# -*- coding: utf8

#########################################
# Modèle d'interface en histogramme
#########################################
import numpy as np
from numpy import linalg as LA

def isPal(nb):
    string = str(nb)
    isp = True
    for y in np.arange(len(string)/2):
        if string[y] != string[-y-1] :
            isp = False
            break
    return isp


NMAX = 1e3
palmax = 0

for i in np.arange(1,NMAX,dtype=int)[::-1] :
    for j in np.arange(1,i)[::-1] :
        if isPal(i*j) and i*j>palmax: 
            palmax = i*j
print palmax
