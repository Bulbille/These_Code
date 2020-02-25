#!/usr/bin/python
# -*- coding: utf8

#########################################
# Modèle d'interface en histogramme
#########################################
import numpy as np
from numpy import linalg as LA

nb = np.loadtxt('nbproj13')
taille = len(nb)

res = 0
for n in nb :
    res += n
print res

