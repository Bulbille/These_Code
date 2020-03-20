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
regex_chain = "mag(\d\.\d+)"
colors=['red','blue','green','rosybrown','sandybrown','indigo','darkmagenta','gold','olivedrab','deepskyblue','lawngreen']
nb_colors = len(colors)

shearspace = []
big_mag = {}
for file in sorted(os.listdir(directory)):
    regex = re.search(regex_chain,file)
    if regex == None :
        continue
    res = re.findall(regex_chain,file)
    shear    = float(res[0])
    if shear not in shearspace:
        shearspace = np.append(shearspace,shear)
    big_mag[shear] = np.loadtxt(directory+file)
shearspace = sorted(shearspace)

def expcor(x,tau,b):
    return np.exp(-x/tau)+b
def expeq(x,tau,a,b):
    return a-b*np.exp(-x/tau)

pmag = plt.subplot(111)
for nt,shear in enumerate(shearspace[:-1]):
    print(nt)
    pmag.plot(big_mag[shear][:,0],big_mag[shear][:,1],color=colors[nt],label="$\\omega="+str(shear)[:4]+'$')

pmag.legend()
pmag.set_xlabel('$y$')
pmag.set_ylabel('$m(y)$')

plt.savefig('profil-mag-ising-shear.pdf')
plt.show()
plt.close()

