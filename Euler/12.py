#!/usr/bin/python
# -*- coding: utf8

#########################################
# Modèle d'interface en histogramme
#########################################
import numpy as np
from numpy import linalg as LA

def nbDiv(x):
    ndiv = 1
    for i in np.concatenate((np.arange(2,x/2+1,dtype=float),[float(x)]),axis=0):
        exp = 0
        while (x/i).is_integer() :
            exp += 1
            x /= i
        ndiv *= exp+1
        if x == 1 :
            break
    return ndiv
    

triangle = 1
increment = 1
MIN = 500

while True :
    increment +=1
    triangle += increment
    res = nbDiv(triangle)
    if res >= MIN :
        print triangle,res
        break
