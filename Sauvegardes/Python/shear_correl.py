#!/usr/bin/python
# -*- coding: utf8

#########################################
# Longueur de corrélation sur les y 
# en fonction de x
# pour chaque température et champ magnétique
#########################################


import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
import scipy as sy
import numpy as np
from pylab import *
import os #pour les chemins de fichiers
from subprocess import call
import sys #pour l'exit et l'argument
import re
import random
r0_255 = lambda: random.randint(0,255)

directory=sys.argv[1]
X=60
Y=200
lx=np.arange(0,X,1)
regex_chain="ttc_([\d.]+)(?:-\d)?(?:-t[\d.]*)?(?:-h-([\d.]+))?"

###############################
# Liste des températures et 
# des champs magnétiques 
###############################
temp = []
champs = []
cisaillement = []

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
big_mag = {}
nb_m = {}

for file in sorted(os.listdir(directory)):
	if re.search(r"correl_y",file) != None :
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

	elif re.search(r"mag_x",file) != None :
		res = re.findall(regex_chain,file)	
		t=float(res[0][0]); h=0
		if bool(res[0][1]):
			h=float(res[0][1])
		data=np.loadtxt(directory+file)
		if (t,h) not in big_mag:
			big_mag[t,h] = data
			nb_m[t,h]=1
		else :
			big_mag[t,h] = np.add(data,big_mag[t,h])
			nb_m[t,h] += 1
		
	else : 
		continue



for points in big_data :
	big_data[points] /= nb[points]
	big_mag[points] /= nb[points]
####################################
# Profils de corrélation 
####################################
def fit_exp(x,l,s,a,h):
	#return -a*np.sin((x-s)/l)+h
	#return -a*np.tanh((x-s)/l)+h
	return a*np.exp(-(x-s)/l)+h

def equation(x,xhi,L):
	def a(xhi,L):
		return (1.-exp(-L/xhi))/(2*sinh(L/xhi))
	def b(xhi,L):
		return (-1.+exp(L/xhi))/(2*sinh(L/xhi))
	return a(xhi,L)*exp(x/xhi)+b(xhi,L)*exp(-x/xhi)

correl = {}

for t in temp:
	for h in champs :
		data=big_data[round(t,3),round(h,3)]
		for x in lx :
			datax=data[:][data[:,0]==x]
			datax=datax[:][datax[:,1]<=30]
			popt=[]
			try:
#				popt,pcov=curve_fit(fit_exp,datax[:Y,0],datax[:Y,1],p0=[4,0,1,0])	
				popt,pcov=curve_fit(fit_exp,datax[:,1],datax[:,2]/np.amax(datax[:,2]),bounds=(0,[9.,0.01,1.1,0.6]))
#				if(x==1 and t==0.8 ) :
#					datax=data[:][data[:,0]==x]
#					plt.plot(datax[:,1],datax[:,2]/np.amax(datax[:,2]))
#					plt.plot(datax[:,1],fit_exp(datax[:,1],popt[0],popt[1],popt[2],popt[3]))
#					plt.plot(datax[:,1],equation(datax[:,1],popt[0],Y))
#					plt.xlabel('Y-distance')
#					plt.ylabel('F_x(y)')
#					plt.title('Correlation function for x='+str(x)+' at t='+str(t)+' and h='+str(h))
#					plt.show()
#					#plt.savefig(directory+'correl_function-t'+str(t)+'h'+str(h)+'x'+str(x)+'.png')
#					plt.close()
			except RuntimeError:
				popt = np.append(popt,Y)
			finally :
				if (t,h,x) not in correl:
					correl[t,h,x] = popt[0]
				else :
					print "++++ Erreur "
					exit()

######################
# Plot
######################

for t in temp :
	for h in champs :
		label=True
		x=[]
		valeur=[]
		for points  in correl:
			if abs(float(points[0])-t)>0.0001 or abs(float(points[1])-h)>0.0001:
			#On n'est pas à la bonne température/champ mag	
				continue
			x=np.append(x,points[2])
			valeur=np.append(valeur,correl[points])

		x,valeur = zip(*sorted(zip(x,valeur)))
		valeur_moy=np.add(valeur[0:X/2],valeur[X/2:X][::-1])
		valeur_moy=tuple([i/2 for i in valeur_moy])
	#	x,valeur = zip(*sorted(zip(x,valeur)))
		plt.plot(x[0:X/2],valeur_moy,label="t="+str(t)+" , h="+str(h))
	plt.legend(loc='upper center',prop={'size':10})
	plt.xlabel('X')
	plt.ylabel('Y-Correlation length')
	plt.title('Y-correlation length')
	plt.savefig(directory+'lc_t'+str(t)+'.pdf')
#	plt.show()
	plt.close()

for h in champs :
	for t in temp :
		label=True
		x=[]
		valeur=[]
		for points  in correl:
			if abs(float(points[0])-t)>0.0001 or abs(float(points[1])-h)>0.0001:
			#On n'est pas à la bonne température/champ mag	
				continue
			x=np.append(x,points[2])
			valeur=np.append(valeur,correl[points])

		x,valeur = zip(*sorted(zip(x,valeur)))
		valeur_moy=np.add(valeur[0:X/2],valeur[X/2:X][::-1])
		valeur_moy=tuple([i/2 for i in valeur_moy])
	#	x,valeur = zip(*sorted(zip(x,valeur)))
		plt.plot(x[0:X/2],valeur_moy,label="t="+str(t)+" , h="+str(h))
	plt.legend(loc='upper center',prop={'size':10})
	plt.xlabel('X')
	plt.ylabel('Y-Correlation length')
	plt.title('Y-correlation length')
	plt.savefig(directory+'lc_h'+str(h)+'.pdf')
#	plt.show()
	plt.close()


####################################
# Fit du profil de magnétisation
####################################
mag = {}
nb_mag = {}

def fit_tanh(x,l,s,a,h):
	return -a*np.tanh((x-s)/l)+h

for t in temp :
	for h in champs :
		data=big_mag[t,h]
		popt=[]
		try:
			popt,pcov=curve_fit(fit_tanh,data[:X/2,0],data[:X/2,1],p0=[3,15,1,0])	
		except RuntimeError :
			popt = np.append(popt,0)
		finally :
	#	print popt, pcov
			if (t,h) not in mag:
				mag[t,h] = popt[0]
			else :
				mag[t,h] += popt[0]

####################################
# Fit du profil en gaussienne
####################################
def fit_gauss(x,l,s,a,h):
	return a*np.exp(-(x-s)**2/(2*l**2))+h

for t in temp :
	champ=[]
	longueurs=[]
	for h in champs :
		label=True
		x=[]
		valeur=[]
		for points  in correl:
			if abs(float(points[0])-t)>0.0001 or abs(float(points[1])-h)>0.0001:
			#On n'est pas à la bonne température/champ mag	
				continue
			x=np.append(x,points[2])
			valeur=np.append(valeur,correl[points])
		x,valeur = zip(*sorted(zip(x,valeur)))
		valeur_moy=np.add(valeur[0:X/2],valeur[X/2:X][::-1])
		valeur_moy=tuple([i/2 for i in valeur_moy])
		popt,pcov=curve_fit(fit_gauss,x[10:20],valeur_moy[10:20],bounds=([0,10,0,0],[8,20,8,0.6]))
		champ=np.append(champ,h)	
		longueurs=np.append(longueurs,np.abs(popt[0]))

	ch_mag=[] ; l_mag=[];
	for points in mag :
		if points[0] != t: 
			continue
		ch_mag=np.append(ch_mag,points[1])
		l_mag=np.append(l_mag,mag[points])
		ch_mag,l_mag = zip(*sorted(zip(ch_mag,l_mag)))
#		if popt[0] <= 0:
#			print popt
#			plt.plot(x[:30],valeur[:30])
#			plt.plot(x[:30],fit_gauss(x[:30],popt[0],popt[1],popt[2],popt[3]))
#			plt.show()
	plt.plot(champ,longueurs,'-+',label="t="+str(t)+" from fits ",color=colors[t])
	plt.plot(ch_mag,l_mag,'--',label="t="+str(t)+" from mag profile",color=colors[t])
plt.legend()
plt.xlabel('Magnetic field')
plt.ylabel('X-Correlation length')
plt.title('X-correlation length after a gaussian fit of the Y-correlation length ')
plt.savefig(directory+'longueur_correl.pdf')
plt.show()
