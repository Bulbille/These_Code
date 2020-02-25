#!/usr/bin/python
# -*- coding: utf8

#########################################
# Modèle d'interface en histogramme
#########################################
import numpy as np
from numpy import linalg as LA

nbfile = 0

pyramid = {}

with open('p067_triangle.txt') as fp :
    for i,line in enumerate(fp):
        pyramid[i] =np.fromstring(line,dtype=int,sep=' ')
    nline = i

print pyramid
for line in np.arange(nline)[::-1]:
    for j,k in enumerate(pyramid[line]):
        pyramid[line][j] += max(pyramid[line+1][j],pyramid[line+1][j+1])
print pyramid
print pyramid[0][0]
