#!/usr/bin/python
# -*- coding: utf8

#########################################
# Modèle d'interface en histogramme
#########################################
import numpy as np
from numpy import linalg as LA

def isTriplet(x,y,z):
    return x**2+y**2 == z**2

MAX = 1000

for a in np.arange(1,MAX)[::-1]:
    for b in np.arange(1,MAX-a)[::-1]:
        c = MAX-b-a
        if isTriplet(a,b,c):
            print a*b*c
