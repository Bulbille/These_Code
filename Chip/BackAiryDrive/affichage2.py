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

import numpy as np
import os #pour les chemins de fichiers
import sys #pour l'exit et l'argument
import re
import scipy as sy
import scipy.special as sp
import scipy.integrate as integrate
from scipy.optimize import curve_fit


directory=sys.argv[1]
regex_chain = "[\d]+.+_mu(\d\.\d+)"
colors=['red','blue','green','rosybrown','sandybrown','indigo','darkmagenta','gold','olivedrab'] ; nb_colors = len(colors)

alpha   =sp.ai_zeros(3)[1]
def p(h,zeros,l,mean):
    nominateur = np.power(sp.airy(abs((h-mean)/l)+zeros[0])[0],2)
    denominateur = 2*integrate.quad(lambda x : np.power(sp.airy(abs(x/l)+zeros[0])[0],2),0,np.inf)[0]
    return nominateur/denominateur


###############################
# Mise en mémoire des données #
###############################
muspace = []
tauxspace = []
fspace = []
big_histo = {}
big_snap = {}
for file in sorted(os.listdir(directory)):
    regex = re.search(regex_chain,file)
    if regex == None : 
            continue
    res = re.findall(regex_chain,file)
    if float(res[0]) not in muspace:
        muspace = np.append(muspace,float(res[0]))
    if re.search('histo',file) :
        big_histo[float(res[0])] = np.loadtxt(directory+file)
    elif re.search('snap',file) :
        big_snap[float(res[0])] = np.loadtxt(directory+file)
muspace = sorted(muspace)
histogramme = plt.subplot(111)

x= np.arange(300,500)
for nm,mu in enumerate(muspace):
    print(mu)
    integraleA = np.sum(big_histo[mu][:,1])
    histogramme.plot(big_histo[mu][:,0],big_histo[mu][:,1]/integraleA,color=colors[nm%nb_colors],label=str(mu))
    popt,pcov =  curve_fit(lambda x,l,mean : p(x,alpha,l,mean), big_histo[mu][:,0],big_histo[mu][:,1]/integraleA,p0=[10,200])
    print(popt)
    histogramme.plot(big_histo[mu][:,0],p(big_histo[mu][:,0],alpha,*popt),'+',color=colors[nm%nb_colors])

histogramme.legend()
histogramme.set_xlim([100,300])

plt.show()
plt.close()


