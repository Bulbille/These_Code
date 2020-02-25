#!/user/bin/python
# -*- coding: utf8

###############################################################
# Affiche l'évolution des profils en fonction de la temperature
# et du champ magnétique/cisaillement/autre paramètre
###############################################################

import matplotlib.pyplot as plt
import scipy as sy
import numpy as np
import os #pour les chemins de fichiers
import sys #pour l'exit et l'argument
import fnmatch
import re

directory='./champ_mag/'

temp = {}
shear = {}
nbfichier = {}
ex = {}
ey = {}
mag = {}

nb_fichier=0

### Magnetisation
for file in sorted(os.listdir(directory)):
        if re.search(".+mag_x",file) == None:
                continue
        data=np.loadtxt(directory+file)
	t = 0
        w = 0
	if re.findall(r'ttc_(\d\.\d*)',file):
		t =  re.findall(r'ttc_(\d\.\d*)',file)[0]
	if re.findall(r'h-(\d\.\d*)-mag_x',file):
		w =  re.findall(r'h-(\d\.\d*)-mag_x',file)[0]

        if t not in temp :
                temp[t]=0
        if w not in shear :
                shear[w]=0

        if (t,w) in mag :
                mag[t,w] += data[:,1]
        else :
                mag[t,w] = data[:,1]
	
## Énergies
for file in sorted(os.listdir(directory)): 
	if re.search(".+energie",file) == None: 
		continue
	data=np.loadtxt(directory+file)
	w = 0
	t = 0
	if re.findall(r'ttc_(\d\.\d*)',file):
		t =  re.findall(r'ttc_(\d\.\d*)',file)[0]
	if re.findall(r'h-(\d\.\d*)-energie',file):
		w =  re.findall(r'h-(\d\.\d*)-energie',file)[0]
	if t not in temp :
		temp[t]=0
	if w not in shear :
		shear[w]=0

	if (t,w) in ex :
		ex[t,w] += data[:,1]
	else : 
		ex[t,w] = data[:,1]

	if (t,w) in ey :
		ey[t,w] += data[:,2]
	else : 
		ey[t,w] = data[:,2]

	if (t,w) in nbfichier :
		nbfichier[t,w] += 1
	else : 
		nbfichier[t,w] = 1

temp = sorted(temp)
shear = sorted(shear)
for t in temp :
	for w in shear :
		ex[t,w] /= nbfichier[t,w]
		ey[t,w] /= nbfichier[t,w]
		mag[t,w] /= nbfichier[t,w]


#########################
#### Données article ####
#########################
#abraham_ex = np.loadtxt('abraham-ex',delimiter=',')
#w1y = np.loadtxt('ey-w01.0',delimiter=',')
#w05y = np.loadtxt('ey-w00.5',delimiter=',')
#abraham_ey = np.loadtxt('abraham-ey',delimiter=',')
#
for t in temp :
	plt.figure(t)
	plt.title("Temperature : " + str(t))
	for w in shear :
		if float(w) > 1 :
			continue
		plt.subplot(311)
		plt.plot(np.linspace(1/len(ex[t,w]),1,len(ex[t,w])),ex[t,w],'-',label='h=' + str(w))
		plt.subplot(312)
		plt.plot(np.linspace(1/len(ey[t,w]),1,len(ex[t,w])),ey[t,w],'-',label='h=' + str(w))
		plt.subplot(313)
		plt.plot(np.linspace(0,len(mag[t,w]),len(mag[t,w])),mag[t,w],'-',label='h=' + str(w))

	plt.subplot(311)
	plt.title("Energy bond x")
	plt.legend(loc=0,prop={'size':10},numpoints=1,ncol=2,labelspacing=0.1,columnspacing=0.5)
	plt.subplot(312)
	plt.title("Energy bond y")
	plt.legend(loc=0,prop={'size':10},numpoints=1,ncol=2,labelspacing=0.1,columnspacing=0.5)
	plt.subplot(313)
	plt.title("Magnetization x")
	plt.legend(loc=0,prop={'size':10},numpoints=1,ncol=2,labelspacing=0.1,columnspacing=0.5)

plt.show()

