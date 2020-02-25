#!/usr/bin/python
# -*- coding: utf8

#########################################
# Modèle d'interface en histogramme
#########################################
import matplotlib.pyplot as plt
fontsize=8
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
LX = 64

directory=sys.argv[1]
regex_chain=".+X([\d]+)F([\d.]+)"
colors=['red','blue','green','rosybrown','sandybrown','indigo','darkmagenta','gold','olivedrab'] ; nb_colors = len(colors)

###############################
# Mise en mémoire des données #
###############################
F = 0
Lspace = []
big_data = {}
for file in sorted(os.listdir(directory)):
    regex = re.search(regex_chain,file)
    if regex == None : 
            continue
    res = re.findall(regex_chain,file)
    if float(res[0][0]) not in Lspace:
        Lspace = np.append(Lspace,float(res[0][0]))
        data=np.loadtxt(directory+file)
        big_data[float(res[0][0])] = data
    F = float(res[0][1])

Lspace  = sorted(Lspace)

############################
# Calcul de l'énergie libre
############################
def linear(x,a,b):
    return a*x+b
def fit_pow(x,a,b,p):
    return a*np.exp(x/p)+b

drive = plt.subplot(211)
energie = plt.subplot(212)

for nb,LX in enumerate(Lspace):
    print nb,LX
    drive.plot(big_data[LX][:,0],big_data[LX][:,1],'+',color=colors[nb],label=str(LX))
    energie.plot(big_data[LX][:,0],big_data[LX][:,3],'+',color=colors[nb],label=str(LX))

    popt,pcov = curve_fit(fit_pow,big_data[LX][:,0],big_data[LX][:,1],p0=[1,18,1])
    print '1  ',popt
    drive.plot(big_data[LX][:,0],fit_pow(big_data[LX][:,0],*popt),'-',color=colors[nb])

    popt,pcov = curve_fit(fit_pow,big_data[LX][:,0],big_data[LX][:,3],p0=[1,1,1])
    print '2  ',popt
    energie.plot(big_data[LX][:,0],fit_pow(big_data[LX][:,0],*popt),'-',color=colors[nb])


drive.set_xlabel('Drive force')
drive.set_ylabel('$\\bar{H}$')
drive.legend(loc='upper left')
energie.set_xlabel('Drive force')
energie.set_ylabel('$E=\sum |h_i-h_{i+1}|$')
energie.legend(loc='upper left')


plt.savefig(directory+'drivemixeddynamics.pdf')
plt.show()
plt.close()


