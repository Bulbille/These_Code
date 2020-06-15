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


def pair(h,zeros,l,mean):
    nominateur = np.power(sp.airy(abs((h-mean)/l)+zeros)[0],2)
    denominateur = 2*integrate.quad(lambda x : np.power(sp.airy(abs(x/l)+zeros)[0],2),0,np.inf)[0]
    return nominateur/denominateur
def impair(h,zeros,l,mean):
    nominateur = (1*(h>0)-1*(h<0))*np.power(sp.airy(abs((h-mean)/l)+zeros)[0],2)
    denominateur = 2*integrate.quad(lambda x : np.power(sp.airy(abs(x/l)+zeros)[0],2),0,np.inf)[0]
    return nominateur/denominateur


N = 5
alpha   =sp.ai_zeros(N)
m = 10

etats    = plt.subplot(111)
hspace = np.linspace(-m,m,200)
for i in np.arange(N):
    if i%2 == 0:
        etats.plot(hspace,pair(hspace,alpha[1][i],1,0),label="$n="+str(i)+"$")
    else:
        etats.plot(hspace,impair(hspace,alpha[0][i],1,0),label="$n="+str(i)+"$")

def gauss(h) :
    return np.exp(-h*h/2)/np.power(2*np.pi,0.5)


etats.plot(hspace,gauss(hspace),color="black",label="Gaussienne")

etats.legend()
etats.set_xlabel('$z$')
etats.set_ylabel('$\psi_n(z)$')
plt.savefig('etats-laser.pdf')
plt.show()
plt.close()


