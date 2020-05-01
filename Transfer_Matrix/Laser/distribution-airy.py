#!/usr/bin/python
# -*- coding: utf8
import matplotlib.pyplot as plt
fsize=10
plt.rc("font",size=fsize)
plt.rc("font",family='serif')

import numpy as np
from scipy.optimize import curve_fit
from numpy import linalg as LA
import scipy as sy
import scipy.special as sp
import scipy.integrate as integrate
import sys
colors=['red','blue','green','rosybrown','sandybrown','gold','olivedrab']
ncol = len(colors)
################################## 
### Déclaration fonctions ########
##################################
def linear(x,a,b):
    return a*x+b
### Définition des matrices de transfert des différents modèles
def Trans(kbt,champ,inter,L):
    d=np.zeros((L,L))
    for y1,i in enumerate(d):
        for y2,j in enumerate(i):
            d[y1][y2] = np.exp(-kbt*( champ*(abs(y1-L/2)+abs(y2-L/2))/2 +inter*abs(y1-y2)) )
    return d
def TransGauss(kbt,champ,inter,L):
    d=np.zeros((L,L))
    for y1,i in enumerate(d):
        for y2,j in enumerate(i):
            d[y1][y2] = np.exp(-kbt*( champ*(abs(y1-L/2)+abs(y2-L/2))/2 +inter*(y1-y2)**2) )
    return d
### Fonction qui calcule directement l'énergie libre
def distrib(ly,beta,h,j,mod):
    if mod == 'SOS' :
        tm = Trans(beta,h,j,ly)
    else :
        tm = TransGauss(beta,h,j,ly)
    w,v = LA.eigh(tm) 
    Z = np.sum(np.power(w,ly))
    ph = np.zeros(ly)
    for l in np.arange(ly):
        ph[l] = np.sum(np.power(w,ly)*np.power(v[l,:],2))
    return ph/Z,-1/np.log(w[-2]/w[-1])

alpha   =sp.ai_zeros(3)[1][0]
alphap   =sp.ai_zeros(3)[0][0]

def p(h,zeros,mean,l):
    nominateur = np.power(sp.airy(abs((h-mean)/l)+zeros)[0],2)
    denominateur = 2*integrate.quad(lambda x : np.power(sp.airy(x/l+zeros)[0],2),0,np.inf)[0]
    return nominateur/denominateur

def pow(x,b,a):
    return b*np.power(x,a)

############################
# Déclaration matrices
J = 1
lMax = 200
lSpace = np.arange(lMax)
fr = 30
xfit = lSpace[int(lMax/2)-fr:int(lMax/2)+fr+1]


fm = ['+','.','v','o','x','^','>','<',0,1,2,3,4,5,6,7,8,9,10,11,'X','1','2','3','4','8']
nfm = len(fm)
tSpace = np.linspace(1,10,5)
bSpace = np.logspace(-3,0,5)

figperp = plt.figure()
perp = figperp.add_subplot(1,1,1)
figpara = plt.figure()
para = figpara.add_subplot(1,1,1)

for i,t in enumerate(tSpace):
    for j,B in enumerate(bSpace):
        print(t,B)
        res = distrib(lMax,1/t,B,J,'Gauss')

        ph = res[0]
        phfit = ph[int(lMax/2)-fr:int(lMax/2)+fr+1]
        popt,pcov= curve_fit(lambda lSpace,l : np.log(p(lSpace,alpha,lMax/2,l)), xfit,np.log(phfit),p0=[3])
        print(popt)
        xiperp = abs(popt[0])
        Sigma = t**2/(2*xiperp**3*B)
        perp.scatter(xiperp,Sigma,marker=fm[j],color=colors[i%ncol],label='$T='+str(t)[:3]+',B='+str(B)[:5]+'$' if j == i else '')

        xipara = res[1]
        Sigma = -B**2/t*xipara**3*(alphap-alpha)**3/2
        para.scatter(xipara,Sigma,marker=fm[j],color=colors[i%ncol],label='$T='+str(t)[:3]+',B='+str(B)[:5]+'$' if j == i else '')
perp.set_ylabel('$\sigma = (2 \\beta^2 B \\xi_\perp^3)^{-1}$')
perp.set_xlabel('$\\xi_\perp$')
para.set_xlabel('$\\xi_\parallel$')
para.set_ylabel('$\sigma = \\frac{\\beta B^2}{2} \left( \\frac{\\alpha_1-\\alpha_0}{\ln(\\frac{\lambda_1}{\lambda_0})} \\right)^3$')
perp.legend()
para.legend()
figperp.savefig('sigma-perp-sos.pdf')
#figperp.show()
figpara.savefig('sigma-para-sos.pdf')
plt.show()

