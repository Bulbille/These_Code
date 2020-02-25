#!/usr/bin/python
# -*- coding: utf8

#########################################
# Modèle d'interface en histogramme
#########################################
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
import scipy as sy
import scipy.special as sp
import scipy.integrate as integrate
import numpy as np
from pylab import *
import os #pour les chemins de fichiers
from subprocess import call
import sys #pour l'exit et l'argument
import re
import random
import scipy.stats as stats
r0_255 = lambda: random.randint(0,255)

directory=sys.argv[1]
regex_chain="lan([\d.]+)(?:-\d)?(?:-t[\d.]*)?(?:-g-([\d.]+))?"

###############################
# Liste des températures et 
# des champs magnétiques 
###############################
temp = []
champs = []

for file in sorted(os.listdir(directory)):
	a=False
	regex = re.search(regex_chain,file)
	if regex == None : 
		continue
	res = re.findall(regex_chain,file)
	if float(res[0][0]) not in temp:
		temp = np.append(temp,float(res[0][0]))
	if  bool(res[0][1]) :
		if float(res[0][1]) not in champs :
			champs = np.append(champs,float(res[0][1]))
	else :
		a=True;
		for i in champs :
			if i<= 0.001:
				a=False
				break
	if a :
		champs = np.append(champs,0.)
colors= {}
for t in temp:
	colors[t] = ('#%02X%02X%02X' % (r0_255(),r0_255(),r0_255()))

############################
# Moyenne des données
############################
big_data = {}
nb = {}

for file in sorted(os.listdir(directory)):
	if re.search(r"histo",file) != None :
		res = re.findall(regex_chain,file)	
		t=float(res[0][0]); h=0
		if bool(res[0][1]):
			h=float(res[0][1])
		data=np.loadtxt(directory+file)
		data=data[:,:][~np.isnan(data[:,2])]
		if (t,h) not in big_data:
			big_data[t,h] = data
			nb[t,h]=1
		else :
			big_data[t,h] = np.add(data,big_data[t,h])
			nb[t,h] += 1
	else : 
		continue


for points in big_data :
	big_data[points] /= nb[points]

###############################
# Calcul des moments
###############################
zeros=sp.ai_zeros(3)
xairy=np.linspace(-3,3,50)
def p(h,chi,zeros):
        nominateur = np.power(sp.airy(abs(h)/chi+zeros[0])[0],2)
        denominateur = 2*integrate.quad(lambda x : np.power(sp.airy(abs(x)/chi+zeros[0])[0],2),0,np.inf)[0]
        return nominateur/denominateur
chi=np.power(integrate.quad(lambda x : np.power(sp.airy(abs(x)+zeros[1][0])[0],2),0,np.inf)[0]/integrate.quad(lambda x : np.power(x*sp.airy(abs(x)+zeros[1][0])[0],2),0,np.inf)[0],0.5)
valeur_airy = p(xairy,chi,zeros[1])

def calcul_moment(x,y):
        moment=[]
        x_bin=abs(x[1]-x[0])
        moment=np.append(moment,np.sum(y)*x_bin)
        moment=np.append(moment,0)
        for index, valeur in np.ndenumerate(y):
                index=index[0]
                moment[1]+= x[index]*valeur*x_bin
        moment[1] /= moment[0]

        moment=np.append(moment,0)
        for index, valeur in np.ndenumerate(y):
                index=index[0]
                moment[2]+=np.power(x[index]-moment[1],2)*valeur*x_bin
        moment[2] = np.power(moment[2]/moment[0],0.5)

        moment=np.append(moment,0)
        for index, valeur in np.ndenumerate(y):
                index=index[0]
                moment[3]+=np.power((x[index]-moment[1])/moment[2],3)*valeur*x_bin
        moment[3] /= moment[0]
        return moment

moments = {}
for t in temp :
	for h in champs :
		x=[];valeur=[];
		for point in big_data[t,h]:
			x=np.append(x,point[0]); valeur=np.append(valeur,point[2])
		moments[t,h] = calcul_moment(x,valeur)
		print "+++++++++++",np.sum(valeur)
		print moments[t,h]
#		x=(x-moments[t,h][1])/moments[t,h][2]; valeur = np.multiply(valeur, moments[t,h][2]/moments[t,h][0])
		print np.sum(valeur)*(x[1]-x[0])
		plt.plot(x,valeur,label="t="+str(t)+",h="+str(h),color=colors[t])
plt.plot(xairy,valeur_airy,".", label="airy", color="black")
plt.plot(x,np.exp(-x*x/2)/np.power(2*np.pi,0.5),"-+",label="gaussienne",color="black")
plt.plot(x,(2/np.pi)/(1+x*x*4),"--",label="lorentz",color="black")
#plt.yscale('log')
plt.legend()
plt.show()
plt.close()
