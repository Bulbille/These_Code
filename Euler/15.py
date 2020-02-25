#!/usr/bin/python
# -*- coding: utf8

#########################################
# Modèle d'interface en histogramme
#########################################
import numpy as np
from numpy import linalg as LA
import math

N = 1000
nb = pow(2,1000)

res = 0
for i,j in enumerate(str(nb)) :
    res += int(j)
print res
