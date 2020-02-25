#!/usr/bin/python
# -*- coding: utf8

############################################
# Grandeurs thermodynamiques à l'équilibre
# Prend les dernières valeurs dans le 
# fichier thermo-kaw-x et en fait une moyenne
############################################

import matplotlib.pyplot as plt
import numpy as np 
import scipy as sy
from pylab import *
import os #pour les chemins de fichiers
import sys
import re

directory='./glau-mag/'

#temp = np.around(np.arange(0.7,1.1,0.01),2)
#temp = np.linspace(0.8,1.05,30)
temp = []
#On va prendre les températures
for file in sorted(os.listdir(directory)):
        if re.search("^thermo_kaw-.$",file) == None:
                continue
	with open(directory+file) as f:
		noms = f.readline().split('\t')
		dtipus = [('beta', sy.float32)] + [('mag', sy.float32)]+ [('energie', sy.float32)]+ [('chi', sy.float32)]+ [('c_v', sy.float32)]
		data = sy.loadtxt(f,delimiter=' ',dtype=dtipus)
	for line in data:
		if line[0] not in temp:
			temp = np.append(temp,line[0])
	
mtot = {}
etot = {}
chitot = {}
cvtot = {}
for t in temp :
	mtot[t] = 0
	etot[t] = 0
	chitot[t] = 0
	cvtot[t] = 0
nb_t={}
for t in temp :
	nb_t[t] = 0


### Thermo

for file in sorted(os.listdir(directory)):
        if re.search("^thermo_kaw-.$",file) == None:
                continue
	data=np.loadtxt(directory+file)
	for t in temp :
		datat=data[:,1:][abs(data[:,0]-t)<=0.001]
		for datat2 in datat:
			datat=datat2[0]
			mtot[t] += datat2[0]
			etot[t]+= datat2[1]
			chitot[t] += datat2[2]
			cvtot[t] += datat2[3]
			nb_t[t]+=1
for t in temp :
	mtot[t] /= nb_t[t]
	etot[t] /= nb_t[t]
	chitot[t] /= nb_t[t]
	cvtot[t] /= nb_t[t]



plt.subplot(221)
plt.plot(mtot.keys(),mtot.values(), "o")
plt.tight_layout()
plt.title("Valeur absolue de la Magnetisation")

plt.subplot(222)
plt.plot(etot.keys(),etot.values(), "o")
plt.tight_layout()
plt.title(u"Énergie interne")

plt.subplot(223)
plt.plot(chitot.keys(),chitot.values(), "o")
plt.tight_layout()
plt.title(u"Susceptibilité magnétique")

plt.subplot(224)
plt.plot(cvtot.keys(),cvtot.values(), "o")
plt.tight_layout()
plt.title(u"Capacité calorifique")

plt.show()
