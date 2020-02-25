#!/usr/bin/python
# -*- coding: utf8

#########################################
# Modèle d'interface en histogramme
#########################################
import matplotlib.pyplot as plt
fontsize=10
plt.rc("font",size=fontsize)
plt.rc("font",family='serif')

import scipy as sy
from scipy.optimize import curve_fit
import numpy as np
from numpy import linalg as LA
import os #pour les chemins de fichiers
import sys #pour l'exit et l'argument
import re

J = 1 ; TC = 2/np.log(1.+np.power(2,0.5)) ; BETA = 1/TC;
LX = 70
Model="C"

directory=sys.argv[1]
regex_chain=Model+".+X[\d.]+Y([\d.]+)"
colors=['red','blue','green','rosybrown','sandybrown','indigo','darkmagenta'] ; nb_colors = len(colors)
Lspace = []

###############################
# Mise en mémoire des données #
###############################
for file in sorted(os.listdir(directory)):
    regex = re.search(regex_chain,file)
    if regex == None : 
            continue
    res = re.findall(regex_chain,file)
    if int(res[0]) not in Lspace:
        Lspace = np.append(Lspace,int(res[0]))

big_data = {}

for file in sorted(os.listdir(directory)):
    if re.search(regex_chain,file) == None :
        continue
    res = re.findall(regex_chain,file)	

    ly=int(res[0])
    try:
        data=np.loadtxt(directory+file)
        big_data[ly] = data
    except: 
        continue
Lspace  = sorted(Lspace)
Hspace  = big_data[Lspace[0]][:,0]
Hmax    = Hspace[-1]

###################################
# Intégration de la magnétisation #
###################################

def intMag(data,hmin,hmax):
    omega = 0
    hn = 0
    for k,h in enumerate(Hspace):
        if h < hmin :
            continue
        hn+=1
        if(h == hmin or h == hmax):
            omega += data[k]
        elif (k%2 == 0 ):
            omega += 2*data[k]
        else:
            omega += 4*data[k]
    omega /= (3*hn)/(hmax-hmin)
    return omega

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

def intDiag(lx,ly,beta,h,j,numero):
    tm = Trans(beta,h,j,ly,numero)
    w,v = LA.eigh(tm) ;
    lZ = 1*np.log(max(w))
    return -1/beta*lZ

def fitexp(x,l,a):
    return a*np.exp(-x/l)
def fitpow(x,l,a):
    return a*np.power(x,-l)

############################
# Calcul de l'énergie libre
############################
for nb,Hmin in enumerate(Hspace):
    if Hmin > 0:
        continue
    FreInf = intMag(big_data[Lspace[-1]][:,1],Hmin,Hmax)
    FreSim = np.empty([np.size(Lspace),2])
    FreTM = np.empty([np.size(Lspace),2])

    for numero,LY in enumerate(Lspace):
        FreSim[numero] = [LY,-intMag(big_data[LY][:,1],Hmin,Hmax)]

        FreTM[numero] = [LY,intDiag(LX,LY,BETA,Hmin,J,Model)]
    print FreSim

    if Model == 'C':
#        FreSim[:,1] -= FreSim[:,0]*Hmax
        FreSim[:,1] = FreSim[:,1]+FreSim[-1,1]
#        FreSim[:,1] = FreSim[:,1]/FreSim[:,0]
#        FreTM[:,1] = FreTM[:,1]-FreTM[-1,1]
#        FreTM[:,1] = FreTM[:,1]/FreTM[:,0]

    plt.plot(FreSim[:-1,0],FreSim[:-1,1],label='$H = '+str(Hmin)+'$',color=colors[nb%nb_colors])
    plt.plot(FreTM[:-1,0],FreTM[:-1,1],'+',color=colors[nb%nb_colors])

#    pope,pcov = curve_fit(fitexp,FreTM[20:-1,0],FreTM[20:-1,1],p0=[4,1])
#    popw,pcov = curve_fit(fitpow,FreTM[:10,0],FreTM[:10,1],p0=[2,1])
#    plt.plot(FreTM[20:-1,0],fitexp(FreTM[20:-1,0],*pope),'o',color=colors[nb%nb_colors])
#    plt.plot(FreTM[:20,0],fitpow(FreTM[:20,0],*popw),'x',color=colors[nb%nb_colors])

#plt.yscale('log')
plt.legend()
plt.savefig(directory+'correlationCh0.3.pdf')
plt.show()
plt.close()


