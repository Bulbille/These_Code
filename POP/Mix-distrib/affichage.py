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
from scipy.optimize import curve_fit


directory=sys.argv[1]
regex_chain = "[\d]+.+_mu(\d\.\d+)_T(\d\.\d+)_F(\d\.\d+)"
colors=['red','blue','green','rosybrown','sandybrown','indigo','darkmagenta','gold','olivedrab'] ; nb_colors = len(colors)

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
    res = re.findall(regex_chain,file)[0]
    if float(res[0]) not in muspace:
        muspace = np.append(muspace,float(res[0]))
    if float(res[1]) not in tauxspace:
        tauxspace = np.append(tauxspace,float(res[1]))
    if float(res[2]) not in fspace:
        fspace = np.append(fspace,float(res[2]))
    
    if re.search('histo',file) :
        big_histo[float(res[0]),float(res[1]),float(res[2])] = np.loadtxt(directory+file)
    elif re.search('snap',file) :
        big_snap[float(res[0]),float(res[1]),float(res[2])] = np.loadtxt(directory+file)
fspace = sorted(fspace)
tauxspace = sorted(tauxspace)
muspace = sorted(muspace)
histogramme = plt.subplot(211)
ecart = plt.subplot(212)

def gauss(x,moy,mean) :
    return 1/(mean*np.power(2*np.pi,0.5))*np.exp(-np.power(x-moy,2)/(2*np.power(mean,2)))

sigmasA  = np.empty([len(fspace),2])
sigmasB = np.empty([len(fspace),2])
sigmasT = np.empty([len(fspace),2])

for nm,mu in enumerate(muspace):
    for nt,taux in enumerate(tauxspace):
        for nf,f in enumerate(fspace):
            print(mu,taux,f)
            courbeX = big_histo[mu,taux,f][:,0]
            courbeA = big_histo[mu,taux,f][:,1]
            courbeB = big_histo[mu,taux,f][:,2]
            courbeT = big_histo[mu,taux,f][:,3]
            integraleA = np.sum(courbeA)
            integraleB = np.sum(courbeB)
            integraleT = np.sum(courbeT)
            histogramme.plot(courbeX,courbeA/integraleA,'-.',color=colors[nf%nb_colors])
            histogramme.plot(courbeX,courbeB/integraleB,'-+',color=colors[nf%nb_colors])
            histogramme.plot(courbeX,courbeT/integraleT,color=colors[nf%nb_colors],label="F "+str(f))

            popt,pcov   = curve_fit(gauss,courbeX,courbeA/integraleA,p0=[20,5])
            sigmasA[nf] = [f,popt[1]]
            popt,pcov   = curve_fit(gauss,courbeX,courbeB/integraleB,p0=[20,5])
            sigmasB[nf] = [f,popt[1]]
            popt,pcov   = curve_fit(gauss,courbeX,courbeT/integraleT,p0=[20,5])
            sigmasT[nf] = [f,popt[1]]

#            if  f == fspace[0] or f == fspace[-1] :
#                snapshot.plot(big_snap[mu,taux,f][:,0],big_snap[mu,taux,f][:,1],'-.',color=colors[nf])
#                snapshot.plot(big_snap[mu,taux,f][:,0],big_snap[mu,taux,f][:,1]+big_snap[mu,taux,f][:,2],label="Mu "+str(mu)+", Taux "+str(taux)+", F "+str(f),color=colors[nf])
#
histogramme.legend()
histogramme.set_xlim([0,75])
histogramme.set_xlabel('$y$')
histogramme.set_ylabel('$p(y)$')

ecart.plot(sigmasA[:,0],sigmasA[:,1],label="Sans drive")
ecart.plot(sigmasB[:,0],sigmasB[:,1],label="Avec drive")
ecart.plot(sigmasT[:,0],sigmasT[:,1],label="Total")
ecart.legend()
ecart.set_xlabel('Drive $f$')
ecart.set_ylabel('Ecart-type du fit en gaussienne')
#snapshot.legend()
#snapshot.set_xlim([0,75])

plt.show()
plt.close()


