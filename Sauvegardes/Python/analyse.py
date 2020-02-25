#!/usr/bin/python
# -*- coding: utf8

#Utiliser l'option "--recalculer" pour forcer à recalculer les fonctions de corrélation
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
import scipy as sy
import numpy as np 
from pylab import *
import os #pour les chemins de fichiers
from subprocess import call
import sys #pour l'exit et l'argument

call(["g++","-o", "Asgard", "correlation.cpp", "../program/hamiltonien.cpp"])
#### FONCTIONS DE FIT #### 
def fit_exp(x,b):
	return np.exp(-x/b)

#### ITÉRATION SUR LES BETA #####

beta = np.array([])
temps_deco = np.array([])
nb_simus = 0


results = open('../simulations/temps_decorrelation','w')

for dossiers in sorted(os.listdir('../simulations')):
	#### ITÉRATION SUR LES TEMPS ####
	if "-" in dossiers and os.path.isdir('../simulations/'+dossiers):
		param,value = dossiers.split("-",1) #value contient la valeur de beta
		if not os.path.isfile('../simulations/'+dossiers+'/correlation') or (len(sys.argv) > 1 and sys.argv[1] == '--recalculer'):
			if nb_simus == 0:

				call(["rm", "../simulations/thermo_glau"])
			#Les dossiers sont formattés "beta-xxx"
			print "............ Beta =" , value
			call(["./Asgard",value]) #Calcul de la fonction de correlation pour chaque Beta
			nb_simus = nb_simus + 1
	else:
		continue

	### LECTURE DU FICHIER DE CORRELATION###
	with open('../simulations/'+dossiers+'/correlation') as f:
		noms = f.readline().split('\t')
		dtipus = [('x', sy.float32)] + [('y', sy.float32)]
		data = sy.loadtxt(f,delimiter=' ',dtype=dtipus)
	x = data['x']
	y = data['y']
	### FIT DE LA FONCTION DE CORRELATION ###
	popt,pcov = curve_fit(fit_exp,x,y,p0=[300])
	beta = np.append(beta,float(value))
	temps_deco = np.append(temps_deco,float(popt[0]))
	results.write(str(value)+' '+str(popt[0])+'\n')

### LECTURE DU FICHIER DES GRANDEURS THERMO###
with open('../simulations/thermo_glau') as f:
	noms = f.readline().split('\t')
	dtipus = [('prout', sy.float32)] + [('mag', sy.float32)]+ [('energie', sy.float32)]+ [('chi', sy.float32)]+ [('c_v', sy.float32)]
	data = sy.loadtxt(f,delimiter=' ',dtype=dtipus)
prout = data['prout']
mag = data['mag']
energie = data['energie']
chi = data['chi']
c_v = data['c_v']

print type(prout[2]), type(beta[2])
prout = 1/prout
beta = 1/beta

plt.subplot(321)
plt.plot(beta,temps_deco)
plt.title(u'Temps de décorrelation')
plt.ylabel(u'Temps')
plt.xlabel(u'Beta')
plt.tight_layout()

plt.subplot(322)
plt.plot(prout,mag)
plt.tight_layout()
plt.title("Valeur absolue de la Magnetisation")

plt.subplot(323)
plt.plot(prout,energie)
plt.tight_layout()
plt.title(u"Énergie interne")

plt.subplot(324)
plt.plot(prout,chi)
plt.tight_layout()
plt.title(u"Susceptibilité magnétique")

plt.subplot(325)
plt.plot(prout,c_v)
plt.tight_layout()
plt.title(u"Capacité calorifique")
plt.show()


results.close()
