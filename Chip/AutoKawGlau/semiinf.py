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

directory=sys.argv[1]
regex_chain=".+X[\d.]+F([\d.]+)"
colors=['red','blue','green','rosybrown','sandybrown','indigo','darkmagenta','gold','olivedrab'] ; nb_colors = len(colors)

###############################
# Mise en mémoire des données #
###############################
Tspace = []
big_data = {}
for file in sorted(os.listdir(directory)):
    regex = re.search(regex_chain,file)
    if regex == None : 
            continue
    res = re.findall(regex_chain,file)
    if float(res[0]) not in Tspace:
        Tspace = np.append(Tspace,float(res[0]))
        data=np.loadtxt(directory+file)
        big_data[float(res[0])] = data

Tspace  = sorted(Tspace)
print Tspace

############################
# Calcul de l'énergie libre
############################
def linear(x,a,b):
    return a*x+b

drive = plt.subplot(311)
energie = plt.subplot(312)
variance = plt.subplot(313)
#coef_var = np.power(5e8/(2e3*LX),0.5)
coef_var = 1e-1

for nb,taux in enumerate(Tspace):
    print nb,taux,np.power(coef_var*taux,0.5)
    #drive.errorbar(big_data[taux][:,0],big_data[taux][:,1],yerr=big_data[taux][:,2]*np.power(coef_var*taux,0.5),fmt='+',color=colors[nb],label='$'+str(100*taux)+'\%$ Glauber')
    drive.plot(big_data[taux][:,0],big_data[taux][:,1],'+',color=colors[nb],label='$'+str(100*taux)+'\%$ Glauber')
    variance.plot(big_data[taux][:,0],big_data[taux][:,2],color=colors[nb],label='$'+str(100*taux)+'\%$ Glauber')
    energie.plot(big_data[taux][:,0],big_data[taux][:,3],'+',color=colors[nb],label='$'+str(100*taux)+'\%$ Glauber')


drive.set_xlabel('Drive force')
drive.set_ylabel('$\\bar{H}$')
drive.legend(loc='upper left')
variance.set_xlabel('Drive force')
variance.legend(loc='upper left')
energie.set_xlabel('Drive force')
energie.set_ylabel('$E=\sum |h_i-h_{i+1}|$')
energie.legend(loc='upper left')


plt.savefig(str(directory[:-1])+'.pdf')
plt.show()
plt.close()


