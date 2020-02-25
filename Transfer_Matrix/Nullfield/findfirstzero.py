#!/usr/bin/python
# -*- coding: utf8
import matplotlib.pyplot as plt
fsize=10
plt.rc("font",size=fsize)
plt.rc("font",family='serif')

import numpy as np
from scipy.optimize import curve_fit
from scipy.optimize import brentq
from numpy import linalg as LA
import sys
colors=['red','blue','green','rosybrown','sandybrown','gold','olivedrab']
nb_col = len(colors)
################################## 
### Déclaration fonctions ########
##################################
### Définition des matrices de transfert des différents modèles

############################
# Déclaration matrices
TC = 2/np.log(1.+np.power(2,0.5)) ; BETA = 1;

#betexp = lambda beta : np.exp(-beta)
#theta = lambda eig,h,r : 1/np.tan(h*eig) - r*np.sin(eig)/(1-r*np.cos(eig))
#eigenvalue = np.linspace(0,0.1,1e4)
#func = theta(eigenvalue,LY,betexp(BETA))
#
#zeros = np.where(abs(func)<1e-3)
#print zeros
#lambdamax = func[zeros[0][0]]
#print eigenvalue[zeros[0]]

def eigenvalue(beta,ly,prec):
    betexp = lambda beta : np.exp(-beta)
    theta = lambda eig,h,r : 1/np.tan(h*eig) - r*np.sin(eig)/(1-r*np.cos(eig))

    X = np.linspace(1e-20,1,pow(10,prec))
    func = theta(X[0],ly,betexp(BETA))
    for i,x in enumerate(X[1:]) :
        functmp = theta(x,ly,betexp(BETA))
        if(func*functmp < 0 ):
            zero = brentq(theta,X[i-1],x,args=(ly,betexp(beta)))
            return -np.log(zero)
        else :
            func = functmp
    return False

LY = np.arange(10,100)
free_energy = np.empty(LY.size)
for i,ly in enumerate(LY) :
    free_energy[i] = eigenvalue(BETA,ly,4)

axfree = plt.subplot(111)
axfree.plot(LY,free_energy)

plt.show()
