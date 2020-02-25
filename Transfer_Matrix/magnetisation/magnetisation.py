#!/usr/bin/python
labelsize=16
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
def calcul_moment(x,y):
    moment=[0,0,0]
    x_bin=1.0
    moment[0] = np.sum(y)*x_bin
    for index, valeur in np.ndenumerate(y):
        index=index[0]
        moment[1]+= x[index]*valeur*x_bin
    moment[1] /= moment[0]
    for index, valeur in np.ndenumerate(y):
        index=index[0]
        moment[2]+=np.power(x[index]-moment[1],2)*valeur
    moment[2] = np.power(moment[2]*x_bin/moment[0],0.5)
    return moment

# Definition de la matrice SOS
def transfert_sos(kbt,champ,inter,L):
    d=np.zeros((2*L+1,2*L+1))
    for y1,i in enumerate(d):
        for y2,j in enumerate(i):
            d[y1][y2] = np.exp(-kbt*champ*(abs(L-y1)+abs(L-y2))/2
                    -kbt*inter*abs(y1-y2) )
    return d

def matrice_phi(L):
    d=np.zeros((2*L+1,2*L+1))
    for y,i in enumerate(d):
        d[y][y] = abs(L-y)
    return d


def integrate_mat(lx,ly,hmin,hmax,hn,beta,j):
    omega = 0
    matphi = matrice_phi(ly)
    for k,h in enumerate(np.linspace(hmin,hmax,hn)):
        tm = transfert_sos(beta,h,j,ly)
        w,v = LA.eigh(tm)
        Z = sum(w**lx)
        tmpow = LA.matrix_power(tm,lx)
        trace = np.trace(np.dot(matphi,tmpow))

        phi = trace/Z
        if(k == 0 or k == hn-1):
            omega += phi
        elif (k%2 == 0 ):
            omega += 2*phi
        else:
            omega += 4*phi
    omega /= (3*hn)/(hmax-hmin)
    return omega

def expo(x,a,n,b):
    return a*x**n+b

############################
# Déclaration matrices
J = 1; BETA = 1
Hn = 100; Hmax = 5; Hmin = 0;

lyx = np.arange(10,100,10).astype(int)
lyx = [50]

for LY in [60] :

    Free_phi = np.zeros(np.size(lyx))
    casimir  = np.empty([np.size(lyx),2])
    for l,LX in enumerate(lyx) :
        Free_phi[l] = integrate_mat(LX,LY,Hmin,Hmax,Hn,BETA,J);
        casimir[l] = [LX,Free_phi[l]]

    plt.plot(casimir[:,0],casimir[:,1],'+',label="LY="+str(LY))

#    popt,pcov = curve_fit(expo,casimir[:,0],casimir[:,1],bounds=([0,-3,0],[10,3,2]))
#    print popt
#    plt.plot(casimir[:,0],expo(casimir[:,0],popt[0],popt[1],popt[2]),label="Fit n="+str(popt[1]))

plt.legend()
plt.ylabel('$F(L)$')
plt.xlabel('Length $LX$')
plt.tight_layout()
plt.savefig('sos.pdf')
plt.show()
plt.close()


