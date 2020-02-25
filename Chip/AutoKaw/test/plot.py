#!/usr/bin/python
# -*- coding: utf8

#########################################
# Modèle d'interface en histogramme
#########################################
import matplotlib.pyplot as plt
fontsize=10
plt.rc("font",size=fontsize)
plt.rc("font",family='serif')
import numpy as np
import sys

data = np.loadtxt(sys.argv[1])
plt.plot(data[:,0],data[:,1])
plt.show()
plt.close()


