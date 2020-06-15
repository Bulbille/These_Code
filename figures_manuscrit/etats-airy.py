#!/usr/bin/python
# -*- coding: utf8

#########################################
# Modèle d'interface en histogramme
#########################################
import matplotlib.pyplot as plt
fontsize=14
labelsize=16
plt.rc("font",size=fontsize)
plt.rc("font",family='serif')

import numpy as np
import re
import scipy as sy
import scipy.special as sp
import scipy.integrate as integrate


def p(h,zeros):
    nominateur = np.power(sp.airy(h+zeros)[0],2)
    denominateur = integrate.quad(lambda x : np.power(sp.airy(x+zeros)[0],2),0,np.inf)[0]
    return nominateur/denominateur

        

N = 5
alpha   =sp.ai_zeros(N)

etats    = plt.subplot(111)

hspace = np.linspace(0,5,400)

for i in np.arange(N):
    etats.plot(hspace,p(hspace,alpha[0][i]),label="$n="+str(i)+"$")

etats.legend(loc='right')
etats.set_xlabel('$z$')
etats.set_ylabel('$\psi_n^2(z)$')
plt.savefig('etats-laser.pdf')
plt.show()
plt.close()


