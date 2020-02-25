#!/usr/bin/python
# -*- coding: utf8

#########################################
# Modèle d'interface en histogramme
#########################################
import numpy as np
from numpy import linalg as LA

nb = np.loadtxt('nbproj8',dtype=str)
nb = str(nb)


LARG = 13

maxi = 0
for i in np.arange(len(nb)-LARG):
    tmp = int(nb[i])
    if tmp == 0 :
        continue
    for j in np.arange(1,LARG):
        if int(nb[i+j]) == 0:
            break
        tmp *= int(nb[i+j])
    if tmp > maxi :
        maxi = tmp
print maxi
