#!/usr/bin/python
# -*- coding: utf8
import matplotlib.pyplot as plt
fontsize=10
plt.rc("font",size=fontsize)
plt.rc("font",family='serif')

import scipy as sy
from scipy.optimize import curve_fit
import numpy as np
from numpy import linalg as LA
import time
import sys
colors=['red','blue','green','rosybrown','sandybrown','gold','olivedrab','palegreen','lightseagreen','darkcyan','royalblue','navy','plum','coral','orange','olive','g','c','dodgerblue','violet','fuchsia','crimson','peru','goldenrod','forestgreen','aquamarine','aqua','blueviolet','pink','tomato','linen','antiquewhite','yellow','greenyellow','lime','cyan','indigo','darkmagenta']


################################## 
### Déclaration fonctions ########
##################################

# Definition de la matrice SOS

def transfert_sos(kbt,champ,inter,L):
    d=np.zeros((2*L+1,2*L+1))
    for y1,i in enumerate(d):
        for y2,j in enumerate(i):
            d[y1][y2] = np.exp(-kbt*champ*(abs(L-y1)+abs(L-y2))/2
                        -kbt*inter*abs(y1-y2) )
    return d

# Definition de la matrice Gaussienne
def matrice_phi(L):
    d=np.zeros((2*L+1,2*L+1))
    for y,i in enumerate(d):
        d[y][y] = abs(L-y)
    return d

def simpson(ajout,i,max):
    if(i == 0 or i == max-1):
        return ajout
    elif (i%2 == 0 ):
        return 2*ajout
    else:
        return 4*ajout

############################
# Déclaration matrices

BETA = 1 ; J = 1; 
LX = 40
LY = 20
h = 0.5

matphi = matrice_phi(2*LY+1)
for l,LX in enumerate([40,100,400,500,700]):
    print l
    tm = transfert_sos(BETA,h,J,2*LY+1)
    w,v = LA.eigh(tm)
    Z = sum(w**LX)
    tmpow = LA.matrix_power(tm,LX)
    trace = np.trace(np.dot(matphi,tmpow))
    plt.scatter(LX,trace/Z)
    print trace/Z

plt.legend()
plt.tight_layout()
plt.show()
plt.close()
