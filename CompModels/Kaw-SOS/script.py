#!/usr/bin/python
# -*- coding: utf8
import matplotlib.pyplot as plt
fontsize=12
plt.rc("font",size=fontsize)
plt.rc("font",family='serif')

import scipy as sy
from scipy.optimize import curve_fit
import numpy as np
import sys
from numpy import linalg as LA
colors=['red','blue','green','rosybrown','sandybrown','gold','olivedrab','palegreen','lightseagreen','darkcyan','royalblue','navy','plum','coral','orange','olive','g','c','dodgerblue','violet','fuchsia','crimson','peru','goldenrod','forestgreen','aquamarine','aqua','blueviolet','pink','tomato','linen','antiquewhite','yellow','greenyellow','lime','cyan','indigo','darkmagenta']


################################## 
### Déclaration fonctions ########
##################################

### Définition des matrices de transfert des différents modèles
def TransA(kbt,champ,inter,L):
    d=np.zeros((L,L))
    for y1,i in enumerate(d):
        for y2,j in enumerate(i):
            d[y1][y2] = np.exp(-kbt*( champ*(y1+y2)/2 +inter*abs(y1-y2)) )
    return d
def TransC(kbt,champ,inter,L):
    d=np.zeros((L,L))
    for y1,i in enumerate(d):
        for y2,j in enumerate(i):
            d[y1][y2] = np.exp(-kbt*(champ*(abs(y1-L/2)+abs(y2-L/2))/2 +inter*abs(y1-y2)) )
    return d
### Définition des matrices de magnetisation des différents modèles
def matA(L):
    d=np.zeros((L,L))
    for y,i in enumerate(d):
        d[y][y] = y
    return d
def matC(L):
    d=np.zeros((L,L))
    for y,i in enumerate(d):
        d[y][y] = abs(L/2-y)
    return d
### Fonction qui intègre de 0 à Bmax.
# numero : 0 pour modèle A, 1 pour modèle B, 2 pour modèle C
### Fonction qui calcule directement l'énergie libre
def intDiag(ly,beta,h,j,numero):
    if numero == 0:
        tm = TransA(beta,h,j,ly)
        matphi = matA(ly)
    else :
        tm = TransC(beta,h,j,ly)
        matphi = matC(ly)
    w,v = LA.eigh(tm) 
    phi = np.dot(np.dot(v[:,-1],matphi),v[:,-1]) 
    return -np.log(w[-1])/beta,phi
############################
# Déclaration matrices
BETA = 1; J = 1; 
LY = 20
mod = 1 #0 = \sum_i h_i, 1 = \sum_i |h_i-L/2|

Hspace = 50
Hmax = 6
Hmin = 0
mags = np.linspace(Hmin,Hmax,Hspace)

Fre = np.empty(Hspace)
Mag = np.empty(Hspace)

for nb,H in enumerate(mags):
    res = intDiag(LY,BETA,-H,J,mod)
    Fre[nb] = res[0]
    Mag[nb] = res[1]
    print(H,res[1])

plt.plot(mags,Fre+Hmax*LY/2,label="Transfer Matrix free energy $F(B)+B L/2$")


Int = np.empty(int(len(mags)/2))
for nb,Hm in enumerate(mags):
    if nb == 0 or nb % 2 != 0 :
        continue
    tm = 0
    hn = len(mags[nb:])
    for kp,h in enumerate(mags[nb:]) :
        k = kp+nb
        if(k == 0 or k == hn-1):
            tm += Mag[k]
        elif (k%2 == 0 ):
            tm += 2*Mag[k]
        else:
            tm += 4*Mag[k]
    tm /= (3*hn)/(Hmax-Hm)
    Int[int(nb/2)] = tm
plt.plot(mags[2::2],Int[1:],'+',label="Transfer Matrix integration, $B_2=6$")

##Integration donnees
data = np.loadtxt(sys.argv[1])
mags = data[:,0]
Int = np.empty(int(len(mags)/2))
for nb,Hm in enumerate(mags):
    if nb == 0 or nb % 2 != 0 :
        continue
    tm = 0
    hn = len(mags[nb:])
    for kp,h in enumerate(mags[nb:]) :
        k = kp+nb
        if(k == 0 or k == hn-1):
            tm += data[k,1]
        elif (k%2 == 0 ):
            tm += 2*data[k,1]
        else:
            tm += 4*data[k,1]
    tm /= (3*hn)/(Hmax-Hm)
    Int[int(nb/2)] = tm
plt.plot(mags[2::2],Int[1:],label="Kawasaki, $B_2=6$")

#data = np.loadtxt(sys.argv[2])
#mags = data[:,0]
#Int = np.empty(int(len(mags)/2))
#for nb,Hm in enumerate(mags):
#    if nb == 0 or nb % 2 != 0 :
#        continue
#    tm = 0
#    hn = len(mags[nb:])
#    for kp,h in enumerate(mags[nb:]) :
#        k = kp+nb
#        if(k == 0 or k == hn-1):
#            tm += data[k,1]
#        elif (k%2 == 0 ):
#            tm += 2*data[k,1]
#        else:
#            tm += 4*data[k,1]
#    tm /= (3*hn)/(Hmax-Hm)
#    Int[int(nb/2)] = tm
#plt.plot(mags[2::2],-Int[1:],label="Glauber, $B_2=6$")


#plt.plot(mags,Mag,'+')
#plt.plot(data[:,0],data[:,1],label="MC")


plt.legend(loc='lower left')
plt.ylabel('$F(B_1)-F(\infty)$')
plt.xlabel('$B_2$')
plt.savefig('integration-free-ene.pdf')
plt.show()
plt.close()


