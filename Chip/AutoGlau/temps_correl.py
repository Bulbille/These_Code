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
import os #pour les chemins de fichiers
import sys #pour l'exit et l'argument
import re

directory=sys.argv[1]
Model="C"
regex_chain=Model+".+N([\d.]+)H0\.3[\d.]+"
colors=['red','blue','green','rosybrown','sandybrown','indigo','darkmagenta']

time = []

for file in sorted(os.listdir(directory)):
    regex = re.search(regex_chain,file)
    if regex == None : 
            continue
    res = re.findall(regex_chain,file)
    if int(res[0]) not in time:
        time = np.append(time,int(res[0]))

############################
# Moyenne des données
############################
big_data = {}

for file in sorted(os.listdir(directory)):
    if re.search(regex_chain,file) == None :
        continue
    res = re.findall(regex_chain,file)	

    t=int(res[0])
    try:
        data=np.loadtxt(directory+file)
        big_data[t] = data
    except: 
        continue

###############################
# Accord TM et SOS
###############################
def fit_exp(t,tau):
    return np.exp(-t/tau)

for j,t in enumerate(sorted(time)) :
    plt.plot(big_data[t][:,0],big_data[t][:,1],color=colors[j],label="Average over "+str(t)+" simulations")

rows= np.where(big_data[200][:,0] < 50)
data = big_data[200][rows]
popt,pcov = curve_fit(fit_exp,data[:,0],data[:,1],p0=[1])
plt.plot(big_data[200][:,0],fit_exp(big_data[200][:,0],*popt),label="Fit with $\\tau = "+str(int(popt[0]))+"$",color='black')

plt.xlim(0,100)
plt.ylim(1e-2,1e0)
plt.yscale('log')
plt.legend(loc='lower left')
plt.savefig(directory+'correlationCh0.3.pdf')
plt.show()
plt.close()


