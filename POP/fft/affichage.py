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

directory=sys.argv[1]
regex_chain = "[\d]+.+_mu(\d\.\d+)_F(\d\.\d+)nb(\d\.\d+)"
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

histogramme = plt.subplot(211)
snapshot = plt.subplot(212)

for nm,mu in enumerate(muspace):
    for nt,taux in enumerate(tauxspace):
        for nf,f in enumerate(fspace):
            integrale = np.sum(big_histo[mu,taux,f][:,1])+np.sum(big_histo[mu,taux,f][:,2])
            histogramme.plot(big_histo[mu,taux,f][:,0],big_histo[mu,taux,f][:,1]/integrale,'-.')
            histogramme.plot(big_histo[mu,taux,f][:,0],big_histo[mu,taux,f][:,2]/integrale,'-+')
            histogramme.plot(big_histo[mu,taux,f][:,0],(big_histo[mu,taux,f][:,1]+big_histo[mu,taux,f][:,2])/integrale,label="Mu "+str(mu)+", Taux "+str(taux)+", F "+str(f),color=colors[nf])
            if nm <= 0 and nt <= 0 and  nf <= 0 :
                snapshot.plot(big_snap[mu,taux,f][:,0],big_snap[mu,taux,f][:,1],'-.',color=colors[nf])
                snapshot.plot(big_snap[mu,taux,f][:,0],big_snap[mu,taux,f][:,1]+big_snap[mu,taux,f][:,2],label="Mu "+str(mu)+", Taux "+str(taux)+", F "+str(f),color=colors[nf])

histogramme.legend()
histogramme.set_xlim([0,12])

snapshot.legend()
snapshot.set_xlim([0,256])

plt.show()
plt.close()


