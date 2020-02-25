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
def trans_sos(nb,m) :
    return m/2-nb
    return 1.*nb*m/(m-1)-1.*m/2

def transfert_sos(kbt,champ,inter,l):
    d=np.zeros((l,l))
    for y1,i in enumerate(d):
        for y2,j in enumerate(i):
            d[y1][y2] = np.exp(-kbt*champ*(abs(trans_sos(y1,l))+abs(trans_sos(y2,l)))
                                -kbt*inter*abs(trans_sos(y1,l)-trans_sos(y2,l)))
    return d

def matrice_phi(l):
    d=np.zeros((l,l))
    for y1,i in enumerate(d):
        for y2,j in enumerate(i):
            if y1==y2 :
                d[y1][y2] = abs(trans_sos(y1,l))
    return d

############################
# Déclaration matrices

data=np.loadtxt(sys.argv[1])
plt.plot(data[:,1],data[:,2],'+',label="Glauber")
data=np.loadtxt(sys.argv[2])
plt.plot(data[:,1],data[:,2],'+',label="Kawasaki")

TC  = 2/np.log(1.+np.power(2,0.5)); BETA = 1/TC ; J = 1; 
Hn = 70; Hmax = 0.4; dh = Hmax/Hn
LX = 40
LY = 10

matrix=np.empty([Hn,2])
matphi = matrice_phi(2*LY+1)
for k,h in enumerate(np.linspace(0,Hmax,Hn)):
    tm = transfert_sos(BETA,h,J,2*LY+1)
    w,v = LA.eigh(tm)
    Z = sum(w**LX)
    tmpow = LA.matrix_power(tm,LX)
    trace = np.trace(np.dot(matphi,tmpow))
    matrix[k] = [h,trace/Z]
plt.plot(matrix[:,0],matrix[:,1],label="Transfer Matrix")

plt.ylabel('$F(T,h,L)$')
plt.xlabel('Magnetic field $H$')
plt.legend()
plt.tight_layout()
plt.savefig('kawasaki-glauber.pdf')
plt.show()
plt.close()
