#!/usr/bin/python
# -*- coding: utf8

#########################################
# Modèle d'interface en histogramme
#########################################
import matplotlib.pyplot as plt
fontsize=9
labelsize=16
plt.rc("font",size=fontsize)
plt.rc("font",family='serif')

import numpy as np
import re
import scipy as sy
import scipy.special as sp
import scipy.integrate as integrate



N = 30
alpha   =sp.ai_zeros(N)


etats    = plt.subplot(111)
hspace = np.arange(N)
espace = np.empty(N)
for i in np.arange(N):
    if i%2 == 0:
        espace[i] = alpha[1][i]
    else:
        espace[i] = alpha[0][i]

etats.plot(hspace,espace,'+',label="$\\alpha_k$")
etats.plot(hspace[:-1],espace[1:]-espace[:-1],'+',label="$\\alpha_{k+1}-\\alpha_k$")

etats.legend()
etats.set_xlabel('$n$')
#etats.set_ylabel('')
plt.savefig('energies-laser.pdf')
plt.show()
plt.close()


