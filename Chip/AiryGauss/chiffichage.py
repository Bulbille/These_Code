#!/usr/bin/python
# -*- coding: utf8

#########################################
# Modèle d'interface en histogramme
#########################################
import numpy as np
import os #pour les chemins de fichiers
import sys #pour l'exit et l'argument
import re
import scipy as sy
import scipy.special as sp
import scipy.integrate as integrate
from scipy.optimize import curve_fit


directory=sys.argv[1]
regex_chain = "[\d]+.+Ttc(\d\.\d+)_mu(\d\.\d+)"
colors=['red','blue','green','rosybrown','sandybrown','indigo','darkmagenta','gold','olivedrab','deepskyblue','lawngreen']
nb_colors = len(colors)

alpha   =sp.ai_zeros(3)[1]
def p(h,zeros,l,mean):
    nominateur = np.power(sp.airy(abs((h-mean)/l)+zeros[0])[0],2)
    denominateur = 2*integrate.quad(lambda x : np.power(sp.airy(abs(x/l)+zeros[0])[0],2),0,np.inf)[0]
    return nominateur/denominateur
def powlaw(x,n,a):
    return np.power(x,n)*a


###############################
# Mise en mémoire des données #
###############################
muspace = []
tempspace = []
big_histo = {}
big_snap = {}
for file in sorted(os.listdir(directory)):
    regex = re.search(regex_chain,file)
    if regex == None : 
            continue
    res = re.findall(regex_chain,file)
    temp    = float(res[0][0])
    mu      = float(res[0][1])
    if mu not in muspace:
        muspace = np.append(muspace,mu)
    if temp not in tempspace:
        tempspace = np.append(tempspace,temp)
    if re.search('histo',file) :
        big_histo[temp,mu] = np.loadtxt(directory+file)
    elif re.search('snap',file) :
        big_snap[temp,mu] = np.loadtxt(directory+file)
muspace = sorted(muspace)
muspace = np.asarray(muspace)
mulen = len(muspace)
tempspace = sorted(tempspace)

xfit= np.arange(300,500)

import matplotlib.pyplot as plt
chiplot = plt.subplot(111)

for nt,temp in enumerate(tempspace) :
    chi     = np.empty(len(muspace))
    for nm,mu in enumerate(muspace):
        if mu < 1e-3:
            continue
        x   = big_histo[temp,mu][:,0]
        hi  = big_histo[temp,mu][:,1]
        integraleA = np.sum(hi)
        hi  /= integraleA
        if mu < 1e-2 : 
            start = 10
        else :
            start = 1
        popt,pcov =  curve_fit(lambda xfit,l,mean : p(xfit,alpha,l,mean), x,hi,p0=[start,400])
        chi[nm] = popt[0]
    ide     = np.where((muspace>1e-2)& ( muspace < 0.8))
    popt,pcov = curve_fit(powlaw,muspace[ide],chi[ide],p0=[-0.5,1])
    hamulen = int(mulen/2)
    print(temp,"&",muspace[hamulen],"&",chi[hamulen],"&",popt[0],"&",2*np.power(chi[hamulen]*np.power(temp,-0.70)*np.power(muspace[hamulen],-popt[0]),3),"\\\\")
    print(temp,"&",muspace[mulen-1],"&",chi[mulen-1],"&",popt[0],"&",2*np.power(chi[mulen-1]*np.power(temp,-0.70)*np.power(muspace[mulen-1],-popt[0]),3),"\\\\")

plt.legend()
plt.show()

