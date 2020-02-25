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
import sys
colors=['red','blue','green','rosybrown','sandybrown','gold','olivedrab']
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
def intMag(lx,ly,hmin,hmax,hn,beta,j,numero):
    omega = 0
    matphi = Mat(ly,numero)
    for k,h in enumerate(np.linspace(hmin,hmax,hn)):
        tm = Trans(beta,h,j,ly)
        w,v = LA.eigh(tm)
        # Limite LX → infty
            # m = <L0|D|L0>
        #phi = np.dot( np.dot(v[:,np.argmax(w)],matphi) , v[:,np.argmax(w)])

        #Pas de limite
        Z = sum(w**lx)
        tmpow = LA.matrix_power(tm,LX)
        phi = np.trace(np.dot(matphi,tmpow))/Z

        if(k == 0 or k == hn-1):
            omega += phi
        elif (k%2 == 0 ):
            omega += 2*phi
        else:
            omega += 4*phi
    omega /= (3*hn)/(hmax-hmin)
    return omega
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
for numero,model in enumerate(['A','B','C']) :
    data = np.loadtxt('test/'+model+'-X'+str(LX)+'Y'+str(LY))
#    data = data[:30] #ne pas prendre toutes les valeurs
    #### Calcul de la magnétisation via TM
    magTM = np.empty(np.size(data[:,0]))
    matphi = Mat(LY,model)
    for k,h in enumerate(data[:,0]) :
        tm = Trans(BETA,h,J,LY,model)
        w,v = LA.eigh(tm)
    #    Z = sum(w**LX)
    #    tmpow = LA.matrix_power(tm,LX)
    #    magTM[k] = np.trace(np.dot(matphi,tmpow))/Z
        magTM[k] = np.dot( np.dot(v[:,np.argmax(w)],matphi) , v[:,np.argmax(w)])
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
    FreTM = np.full(np.size(mags),intDiag(LX,LY,BETA,0,J,model))
    for nb,Hmax in enumerate(mags):
        FreTM[nb]  -= intDiag(LX,LY,BETA,Hmax,J,model)

    plt.subplot(211)
    plt.plot(mags,-FreSim,label='Simulation '+model,color=colors[numero])
    plt.plot(mags,FreTM,'+',label='TM '+model,color=colors[numero])
#    plt.plot(mags,-FreMag,'-',label='Int '+model,color=colors[numero])
    plt.subplot(212)
    plt.errorbar(data[:,0],data[:,1],yerr=data[:,2],label='Simulation '+model,color=colors[numero])
    plt.plot(data[:,0],magTM,'+',label='TM '+model,color=colors[numero])

#    if model == 'C':
#        popt,pcov = curve_fit(linear,mags,FreTM,p0=[LY,0])
#        plt.plot(mags,linear(mags,*popt),'x',color=colors[numero],label='Fit $F(B) \propto '+str(round(popt[0],0))+'B$')

### Plot des énergies libre
plt.subplot(211)
plt.xlabel('$B^\\ast$')
plt.ylabel('$F(0)-F(B^{\\ast}) = - \int_0^{B^{\\ast}} m(B) dB$')
#plt.xlim(-0.001,0.3)
plt.legend()
### Plot des magnétisations
plt.subplot(212)
plt.xlim(-0.005,0.3)
plt.xlabel('$B^\\ast$')
plt.ylabel('$m(B)$')
#plt.yscale('log')
plt.legend()

plt.savefig('ModKaw.pdf')
plt.show()


