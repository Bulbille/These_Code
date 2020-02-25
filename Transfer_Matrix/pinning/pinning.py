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
colors=['red','blue','green','rosybrown','sandybrown','gold','olivedrab','palegreen','lightseagreen','darkcyan','royalblue','navy','plum','coral','orange','olive','g','c','dodgerblue','violet','fuchsia','crimson','peru','goldenrod','forestgreen','aquamarine','aqua','blueviolet','pink','tomato','linen','antiquewhite','yellow','greenyellow','lime','cyan','indigo','darkmagenta']


################################## 
### Déclaration fonctions ########
##################################

# Definition de la matrice SOS
def transfert_sos(kbt,champ,inter,L,p):
    d=np.zeros((2*L+1,2*L+1))
    for y1,i in enumerate(d):
        for y2,j in enumerate(i):
            d[y1][y2] = np.exp( -kbt*inter*abs(y1-y2) 
                                +kbt*p*( (y1==0)+(y2==0) )/2 )
    return d

############################
# Déclaration matrices
LX = 50
TC = 2/np.log(1.+np.power(2,0.5));
J =1; 
BETA = 1/TC ; R = np.exp(-BETA)
Pc = -np.log(1-R)
Pn = 4;  
Pecart = 0.5
Pins = np.append(- np.arange(Pn)[::-1][:-1],np.arange(Pn)) / (1.*Pn) * Pecart + Pc

Lmin = 5 ; Lmax = 50 ; Lys = np.arange(Lmin,Lmax) ;  Ln = np.size(Lys)

for nbp,pin in enumerate(Pins):
    print pin
    part = np.empty(Ln)
    for n,LY in enumerate(Lys):
        tm = transfert_sos(BETA,0,J,LY,pin)
        w,v = LA.eigh(tm)
        part[n] = - np.log(max(w))/BETA
    if pin == Pc:
        plt.plot(Lys,part,'o',label="Critical P="+str(pin))
    else :
        plt.plot(Lys,part,label=str(pin))


plt.xlabel('$L_Y$')
plt.ylabel('$F(L) = - \\frac{\ln(\lambda_{max})}{\\beta}$')
plt.legend()
plt.show()
plt.savefig('courbe.pdf')


