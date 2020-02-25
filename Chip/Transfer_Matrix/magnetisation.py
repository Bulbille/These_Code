#!/usr/bin/python
# -*- coding: utf8
import matplotlib.pyplot as plt
fontsize=6
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

def transfert_sos(kbt,champ,inter,L):
    d=np.zeros((2*L+1,2*L+1))
    for y1,i in enumerate(d):
        for y2,j in enumerate(i):
            d[y1][y2] = np.exp(-kbt*champ*(abs(L-y1)+abs(L-y2))/2
                                -kbt*inter*abs(y1-y2) )
    return d

# Definition de la matrice Gaussienne
def transfert_gauss(kbt,champ,inter,L):
    d=np.zeros((2*L+1,2*L+1))
    for y1,i in enumerate(d):
        for y2,j in enumerate(i):
            d[y1][y2] = np.exp(-kbt*champ*(abs(L-y1)+abs(L-y1))
                                -kbt*inter*(y1-y2)**2)
    return d

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
TC  = 2/np.log(1.+np.power(2,0.5)); BETA = 1/TC ; J = 1; 
Hn = 40; Hmax = 0.4; Hmin=0; dh = (Hmax-Hmin)/Hn
LX = 2000

minl = 20; maxl = 40
lys = np.arange(minl,maxl)
lyp = [60,61]


for LX in [40,100,300,600]:
    print LX

    for nbt,t in enumerate(np.linspace(1,3,1)):
        Free_phi = np.zeros(np.size(lys))
        Free_phip = np.zeros(np.size(lyp))
        casimir = np.empty([maxl-minl-1,2])
        BETA = 1./t

        for l,LY in enumerate(lys) :
            matphi = matrice_phi(LY)
            for k,h in enumerate(np.linspace(Hmin,Hmax,Hn)):
                tm = transfert_sos(BETA,h,J,LY)
                tmpow = LA.matrix_power(tm,LX)

                w,v = LA.eigh(tm)
                Z = sum(w**LX)
                phi = np.trace(np.dot(matphi,tmpow)) / Z
#                print LX,LY,Z

                Free_phi[l] += simpson(phi,k,Hn)*dh/3.

        for l,LY in enumerate(lyp) :
            matphi = matrice_phi(2*LY+1)
            for k,h in enumerate(np.linspace(0,Hmax,Hn)):
                tm = transfert_gauss(BETA,h,J,2*LY+1)
                tmpow = LA.matrix_power(tm,LX)

                w,v = LA.eigh(tm)
                Z = sum(w**LX)
                phi = np.trace(np.dot(matphi,tmpow)) / Z

                Free_phip[l] += simpson(phi,k,Hn)*dh/3.
        casimir_p = (Free_phip[0]-Free_phip[1])

        for l,LY in enumerate(lys) :
            if LY < maxl-1 :
                casimir[l] = [LY, - LY**0*( (Free_phi[l]-Free_phi[l+1])- casimir_p )] 
#                casimir[l] = [LY, - (Free_phi[l]-Free_phi[l+1])]
        plt.plot(casimir[:,0],casimir[:,1],'+',label="LX="+str(LX))

plt.legend()
plt.ylabel('$LY^3 (F(L)-F(L+1)$')
plt.xlabel('Length $LY$')
plt.title('No scaling')
plt.tight_layout()
plt.savefig('casimir.pdf')
plt.show()
plt.close()
