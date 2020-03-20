#!/usr/bin/python
# -*- coding: utf8


import numpy as np
import matplotlib.pyplot as plt
fontsize= 12
plt.rc("font",size=fontsize)
plt.rc("font",family='serif')
plt.rc("text",usetex=True) #Latex
import sys
import os
import re
from scipy.optimize import curve_fit

#####################################
# Affiche l'état d'un système de 
# de taille N*M 
#####################################
LY  = 24
J   = 1

directory=sys.argv[1]
regex_chain = "(\d\.\d+)"
colors=['red','blue','green','rosybrown','sandybrown','indigo','darkmagenta','gold','olivedrab','deepskyblue','lawngreen']
nb_colors = len(colors)

tempspace = []
big_eq = {}
big_cor = {}
for file in sorted(os.listdir(directory)):
    regex = re.search(regex_chain,file)
    if regex == None :
        continue
    res = re.findall(regex_chain,file)
    temp    = float(res[0])
    if temp not in tempspace:
        tempspace = np.append(tempspace,temp)
    if re.search('eq',file) :
        big_eq[temp] = np.loadtxt(directory+file)
    elif re.search('cor',file) :
        big_cor[temp] = np.loadtxt(directory+file)
tempspace = sorted(tempspace)

def expcor(x,tau,b):
    return np.exp(-x/tau)+b
def expeq(x,tau,a,b):
    return a-b*np.exp(-x/tau)

pequi = plt.subplot(211)
pcore = plt.subplot(212)
for nt,temp in enumerate(tempspace):
    print(nt,temp)
    popt,pcov = curve_fit(expeq,big_eq[temp][:100,0],big_eq[temp][:100,2],p0=[1,1,10])
    pequi.plot(big_eq[temp][:,0],big_eq[temp][:,2],color=colors[nt],label="$T="+str(temp)[:4]+'$,$\\tau_{eq} = '+str(int(popt[0]))+'$')
    pequi.plot(big_eq[temp][:,0],expeq(big_eq[temp][:,0],*popt),'-.',color=colors[nt])

    popt,pcov = curve_fit(expcor,big_cor[temp][:200,0],big_cor[temp][:200,2],p0=[10,0])
    print(popt)
    pcore.plot(big_cor[temp][:,0],big_cor[temp][:,2],color=colors[nt],label="$T="+str(temp)[:4]+'$,$\\tau_{cor} = '+str(int(popt[0]))+'$')

pequi.set_xlim([0,1000])
pequi.set_xlabel('$t$')
pequi.set_ylabel('$E(t)$')
pequi.legend()

pcore.set_xlim([0,500])
pcore.set_xlabel('$t$')
pcore.set_ylabel('$C(t)/C(0)$')
pcore.legend()

plt.savefig('sos-glau-eq-cor.pdf')
plt.show()
plt.close()

