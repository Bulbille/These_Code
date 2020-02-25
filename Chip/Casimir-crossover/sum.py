#!/usr/bin/python
# -*- coding: utf8

#########################################
# Modèle d'interface en histogramme
#########################################
import matplotlib.pyplot as plt
import scipy as sy
import scipy.special as sp
import scipy.integrate as integrate
import numpy as np
import numpy.linalg as linalg
from pylab import *
import os #pour les chemins de fichiers
import sys #pour l'exit et l'argument
import re

file=sys.argv[1]
data = np.loadtxt(file)
print np.sum(data[:,1])
