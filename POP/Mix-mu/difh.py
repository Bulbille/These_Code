#!/usr/bin/python
# -*- coding: utf8

#########################################
# Modèle d'interface en histogramme
#########################################
import matplotlib.pyplot as plt
fontsize=9
labelsize=16
plt.rc("font",size=fontsize)
plt.rc("font",family='serif')

import scipy as sy
from scipy.optimize import curve_fit
import numpy as np
from numpy import linalg as LA
import math as mt
import os #pour les chemins de fichiers
import sys #pour l'exit et l'argument
import re

TC = 2/np.log(1.+np.power(2,0.5));BETA = 1/TC
LY = 100; J = 1;

directory=sys.argv[1]
taux= 1.00
tauxreg = "%0.2f" % (taux,)
regex_chain = "X[\d]+.+_mu(\d\.\d+)_F"+tauxreg
colors=['red','blue','green','rosybrown','sandybrown','indigo','darkmagenta','gold','olivedrab'] ; nb_colors = len(colors)

def lnfact(n):
    if(n<25):
        return mt.log(mt.factorial(n))
    else:
        return n*mt.log(n)-n

def TransPOP(kbt,inter,L,mu):
    d=np.zeros((L+1,L+1))
    for y1,i in enumerate(d):
        for y2,j in enumerate(i):
            d[y1][y2] = np.exp(-kbt*(-mu*(y1+y2)/2+inter*abs(y1-y2))-(lnfact(y1)+lnfact(y2))/2 )
    return d

### Définition des matrices de magnetisation des différents modèles
def Mat(L):
    d=np.zeros((L+1,L+1))
    for y,i in enumerate(d):
        d[y][y] = y
    return d

### Fonction qui calcule directement l'énergie libre
def intDiag(ly,beta,j,mu):
    matphi = Mat(LY)
    matphi2 = np.square(matphi)
    matphi3 = np.power(matphi,3)
    
    #Dérivée à l'ordre un avec 5 points E = -d(lnZ)/dB
    db = 0.0001
    i=2
    f = np.empty(2*i+1)
    for n in np.arange(2*i+1):
        tm = TransPOP(beta+(n-i)*db,j,ly,mu)
        w,v = LA.eigh(tm) 
        f[n] = mt.log(w[-1])
    E = -(1*f[i-2]-8*f[i-1]+0*f[i+0]+8*f[i+1]-1*f[i+2])/(12*1.0*db**1)/beta

    tm = TransPOP(beta,j,ly,mu)
    w,v = LA.eigh(tm) 

    phi =  np.dot(np.dot(v[:,-1],matphi),v[:,-1])
    phi2 =  np.dot(np.dot(v[:,-1],matphi2),v[:,-1])
    sigma = (phi2-phi**2)**0.5
    gamma = 0
    return [phi,sigma,gamma,E]
###############################
# Mise en mémoire des données #
###############################
muspace = []
big_data = {}
for file in sorted(os.listdir(directory)):
    regex = re.search(regex_chain,file)
    if regex == None : 
            continue
    res = re.findall(regex_chain,file)
    if float(res[0]) not in muspace:
        muspace = np.append(muspace,float(res[0]))
        data=np.loadtxt(directory+file)
        big_data[float(res[0])] = data

muspace  = sorted(muspace)
print muspace


if taux == 0 :
    energie = plt.subplot(311)
    variance = plt.subplot(312)
    skew = plt.subplot(313)
else :
    drive = plt.subplot(221)
    energie = plt.subplot(222)
    variance = plt.subplot(223)
    skew = plt.subplot(224)
derivene = energie.twinx()
mm = 10

for nb,mu in enumerate(muspace):

    print "** Potentiel chimique ",mu
    fspace = big_data[mu][:,0]
    donne  = np.where(fspace<18)
    fspace = fspace[donne]
    totmag = big_data[mu][donne,1][0]
    totvar = big_data[mu][donne,2][0]
    totske = big_data[mu][donne,3][0]
    totene = big_data[mu][donne,4][0]
    # Diagonalisation TM
    res = intDiag(LY,BETA,J,np.linspace(0,40,mu))
    print res
    # Plot moyenne
    if taux > 0 :
        drive.plot(fspace,totmag,color=colors[nb])
    # Plot variance
    variance.plot(fspace,totvar,color=colors[nb])
    # Plot asymétrie
    skew.plot(fspace,totske,color=colors[nb])
    # Plot énergie + dérivée
    energie.plot(fspace,totene,color=colors[nb])
    derivene.plot(fspace,np.gradient(totene/totene[0]),'-.',color=colors[nb])
#        if taux > 0 :
#            popt,pcov = curve_fit(linear,fspace[:mm],totmag[:mm],p0=[1,1,2])
#            drive.plot(fspace[:mm],linear(fspace[:mm],*popt),'+',color=colors[nb])
#            print "Mag", B,popt
#        popt,pcov = curve_fit(linear,fspace[:mm],totvar[:mm]**2/totmag[:mm]**2,p0=[1,1,2])
#        variance.plot(fspace[:mm],linear(fspace[:mm],*popt),'+',color=colors[nb])
#        print "Var", B,popt
#        popt,pcov = curve_fit(linear,fspace[:mm],totene[:mm]/totene[0],p0=[1,1,2])
#        energie.plot(fspace[:mm],linear[:mm](fspace,*popt),'+',color=colors[nb])
#        print "Ene", B, popt
#    except :
#        pass

if taux > 0 :
    drive.set_xlabel('$f$')
    drive.set_ylabel('$\\bar{B}$')
    drive.legend(loc='upper left')
variance.set_xlabel('$f$')
energie.set_xlabel('$f$')
skew.set_xlabel('$f$')

variance.legend(loc='upper left')
energie.legend(loc='upper left')
skew.legend(loc='lower right')

energie.set_ylabel('$E$',fontsize=labelsize)
variance.set_ylabel('$\sigma$',fontsize=labelsize)
skew.set_ylabel('$\gamma$',fontsize=labelsize)

#ligne pointillé pour montrer 2J et 4J
energie.axvline(2*J,linestyle='-.',color='black')
energie.axvline(4*J,linestyle='-.',color='black')
variance.axvline(2*J,linestyle='-.',color='black')
variance.axvline(4*J,linestyle='-.',color='black')
skew.axvline(2*J,linestyle='-.',color='black')
skew.axvline(4*J,linestyle='-.',color='black')
#variance.plot(x1,y1,'-.',x2,y2,'-.',color='black')
#skew.plot(x1,y1,'-.',x2,y2,'-.',color='black')

#variance.set_ylim([2.5,4.5])
#energie.set_ylim([1,8])
#skew.set_ylim([-6,-1])


plt.savefig(str(directory)+'/taux'+str(taux)+'.pdf')
plt.show()
plt.close()


