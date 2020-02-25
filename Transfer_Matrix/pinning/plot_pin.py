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
TC = 2/np.log(1.+np.power(2,0.5));
J =1; 
BETA = 1/TC ; R = np.exp(-BETA)
Pc = -np.log(1-R)
Pn = 7; Pmax = 3; Pmin = 0; Pins = np.linspace(Pmin,Pmax,Pn)

LY = 15
LX = 50

for nbp,pin in enumerate(Pins):
    tm = transfert_sos(BETA,0,J,LY,pin)
    w,v = LA.eigh(tm)
    Z = sum(w**LX)
    h = v[:,np.argmax(w)]**2

    #Distribution
    plt.subplot(211)
    plt.plot(np.arange(0,2*LY+1),h,label="P="+str(pin))
    #h_0
    plt.subplot(212)
    plt.scatter(pin,h[0])

plt.subplot(211)
plt.xlabel('$L_Y$')
plt.ylabel('$p(h) = <\lambda_{max} | h >^2$')
plt.title('Height distribution for $L_Y='+str(LY)+'$ and $L_X='+str(LX)+'$')
plt.legend()

plt.subplot(212)
plt.xlabel('Pinning force')
plt.ylabel('$p(h=0) = <\lambda_{max} | 0 >^2$')
plt.title('Height at the pinning site')
plt.legend()


plt.show()
plt.savefig('height_distrib.pdf')


