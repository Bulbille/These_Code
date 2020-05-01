#!/usr/bin/python
# -*- coding: utf8
import matplotlib.pyplot as plt
fsize=10
plt.rc("font",size=fsize)
plt.rc("font",family='serif')

import scipy as sy
from scipy.optimize import curve_fit
import numpy as np
from numpy import linalg as LA
import sys
colors=['red','blue','green','rosybrown','sandybrown','gold','olivedrab']
fmt=['-','+','o']
################################## 
### Déclaration fonctions ########
##################################
def linear(x,a,b):
    return a*x+b
### Définition des matrices de transfert des différents modèles
def TransA(kbt,champ,inter,L):
    d=np.zeros((2*L+1,2*L+1))
    for y1,i in enumerate(d):
        for y2,j in enumerate(i):
            d[y1][y2] = np.exp(-kbt*champ*((L-y1)+(L-y2))/2 -kbt*inter*abs(y1-y2) )
    return d
def TransB(kbt,champ,inter,L):
    d=np.zeros((2*L+1,2*L+1))
    for y1,i in enumerate(d):
        for y2,j in enumerate(i):
            d[y1][y2] = np.exp(-kbt*( champ*(abs(L-y1)+abs(L-y2))/2 +inter*abs(y1-y2)) )
    return d
def TransC(kbt,champ,inter,L):
    d=np.zeros((2*L+1,2*L+1))
    for y1,i in enumerate(d):
        for y2,j in enumerate(i):
            d[y1][y2] = np.exp(-kbt*( -champ*(abs(L-y1)+abs(L-y2))/2 +inter*abs(y1-y2)) )
    return d
def Trans(kbt,champ,inter,L,numero):
    if numero == 'A':
        return TransA(kbt,champ,inter,L)
    elif numero =='B' :
        return TransB(kbt,champ,inter,L)
    else :
        return TransC(kbt,champ,inter,L)
### Définition des matrices de magnetisation des différents modèles
def matA(L):
    d=np.zeros((2*L+1,2*L+1))
    for y,i in enumerate(d):
        d[y][y] = (L-y)
    return d
def matB(L):
    d=np.zeros((2*L+1,2*L+1))
    for y,i in enumerate(d):
        d[y][y] = abs(L-y)
    return d
def matC(L):
    d=np.zeros((2*L+1,2*L+1))
    for y,i in enumerate(d):
        d[y][y] = -abs(L-y)
    return d
def Mat(L,numero):
    if numero == 'A':
        return matA(L)
    elif numero =='B' :
        return matB(L)
    else:
        return matC(L)
### Fonction qui intègre de 0 à Bmax.
# numero : 0 pour modèle A, 1 pour modèle B, 2 pour modèle C
### Fonction qui calcule directement l'énergie libre
def intDiag(lx,ly,beta,h,j,numero):
    tm = Trans(beta,h,j,ly,numero)
    w,v = LA.eigh(tm) ; 
    lZ = 1*np.log(max(w))
    return -1/beta*lZ
############################
# Déclaration matrices
J = 1 ; TC = 2/np.log(1.+np.power(2,0.5)) ; BETA = 1/TC;
LX = 70
LY = 30


####### Données #####
axmag = plt.subplot(311)
axfree = plt.subplot(313)
axnoise = plt.subplot(312)
axqual = axnoise.twinx()

for numero,model in enumerate(['A','B','C']) :
    try:
        data = np.loadtxt('data/'+model+'-X'+str(LX)+'Y'+str(LY))
        data = data[:50]
        #### Calcul de la magnétisation via TM
        magTM = np.empty(np.size(data[:,0]))
        matphi = Mat(LY,model)
        magTM2 = np.empty(np.size(data[:,0]))
        matphi2 = np.power(Mat(LY,model),2)
        for k,h in enumerate(data[:,0]) :
            tm = Trans(BETA,h,J,LY,model)
            w,v = LA.eigh(tm)
            magTM[k] = np.dot( np.dot(v[:,np.argmax(w)],matphi) , v[:,np.argmax(w)])
            magTM2[k] = pow(np.dot( np.dot(v[:,np.argmax(w)],matphi2) , v[:,np.argmax(w)]) - pow(magTM[k],2) , 0.5)
        ##### Calcul énergie libre via intégration de la magnétisation TM et Simu
        mags = np.empty(np.size(data[:,0])/2)
        FreSim = np.zeros(np.size(mags))
        FreMag = np.zeros(np.size(mags))
        for k in np.arange(np.size(data[:,0])) :
            if k%2 != 0 or k==0:
                continue
            mags[k/2] = data[k,0]
            FreSim[k/2] = FreSim[k/2-1]
            FreSim[k/2] += (mags[k/2]-mags[k/2-1])/6 * ( data[k-2,1] + 4*data[k-1,1] + data[k,1] )

            FreMag[k/2] = FreMag[k/2-1]
            FreMag[k/2] += (mags[k/2]-mags[k/2-1])/6 * ( magTM[k-2] + 4*magTM[k-1] + magTM[k] )
        #### Calcul énergie libre par F(0)-F(B)
        FreTM = np.empty(np.size(mags))
        FreTM[:] = intDiag(LX,LY,BETA,0,J,model)
        for nb,Hmax in enumerate(mags):
            FreTM[nb]  -= intDiag(LX,LY,BETA,Hmax,J,model)
        ### Plot de l'energie libre
        axfree.plot(mags,-FreSim,label='Simulation '+model,color=colors[numero])
        axfree.plot(mags,FreTM,'+',label='TM '+model,color=colors[numero])
        ### Plot de la magnétisation
        axmag.errorbar(data[:,0],data[:,1],yerr=data[:,2],label='Simulation '+model)#,color=colors[numero])
        axmag.errorbar(data[:,0],magTM,yerr=magTM2,marker='o',label='TM '+model)#,color=colors[numero])
        ### Plot du bruit et qualité
        axnoise.plot(data[:,0],abs(magTM2),fmt[numero],label='Noise '+model,color='blue')
        axqual.plot(data[:,0],abs(magTM2/magTM),fmt[numero],label='Quality '+model,color='red')

    except : 
        pass


### Plot des énergies libre
axfree.set_xlabel('$B^\\ast$')
axfree.set_ylabel('$F(0)-F(B^{\\ast}) = - \int_0^{B^{\\ast}} m(B) dB$')
axfree.set_xlim(-0.005,0.3)
axfree.legend(fontsize=fsize-2,loc="center right")
### Plot des magnétisations
axmag.set_xlim(-0.005,0.3)
axmag.set_xlabel('$B^\\ast$')
axmag.set_ylabel('$m(B)$')
axmag.legend(fontsize=fsize-2,loc="center right")

###  Plot du bruit
axnoise.set_xlim(-0.005,0.3)
axnoise.legend(fontsize=fsize-2,loc="upper right")
axnoise.tick_params(axis='y', labelcolor='blue')
### Plot sur le bruit de la qualité 
axqual.legend(fontsize=fsize-2,loc="lower right")
axqual.tick_params(axis='y', labelcolor='red')
#axqual.set_ylim(0,14)

plt.savefig('ModGlau.pdf')
plt.show()


