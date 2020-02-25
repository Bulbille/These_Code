#60!/usr/bin/python
# -*- coding: utf8

#########################################
# Longueur de corrélation sur les y 
# en fonction de x
# pour chaque température et champ magnétique
#########################################


import matplotlib.pyplot as plt
import itertools
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
r0_255 = lambda: random.randint(0,255)
np.seterr(all='warn')

directory=sys.argv[1]
LX=80
LY=50
lx=np.arange(0,LX,1)
regex_chain="ttc_([\d.]+)(?:-\d)?(?:-t[\d.]*)?(?:-h-([\d.]+))?"

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
big_mag = {}
nb_m = {}

for file in sorted(os.listdir(directory)):
	if re.search(r"pdf",file) != None :
		continue
	if re.search(r"correl",file) != None :
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

	elif re.search(r"magy",file) != None :
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
			datax=datax[:][datax[:,1]<=LX/2]
			popt=[]
			try:
#				popt,pcov=curve_fit(fit_exp,datax[:LY,0],datax[:LY,1],p0=[4,0,1,0])	
				popt,pcov=curve_fit(fit_exp,datax[:,1],datax[:,2]/np.amax(datax[:,2]),bounds=([0,0,0.8,-1],[9.,0.1,1.1,1]))
#				if(x==0 or x==LX/4) :
#					datax=data[:][data[:,0]==x]
#					plt.plot(datax[:,1],datax[:,2]/np.amax(datax[:,2]), label="data")
#					plt.plot(datax[:,1],fit_exp(datax[:,1],popt[0],popt[1],popt[2],popt[3]))
#					plt.plot(datax[:,1],equation(datax[:,1],popt[0],LY))
#					plt.xlabel('Y-distance')
#					plt.ylabel('F_x(y)')
#					plt.legend();
#					plt.title('Correlation function for x='+str(x)+' at t='+str(t)+' and h='+str(h))
#					plt.show()
#					#plt.savefig(directory+'correl_function-t'+str(t)+'h'+str(h)+'x'+str(x)+'.png')
#					plt.close()
			except :
				popt=np.append(popt,0)
			finally :
				if (t,h,x) not in correl:
					correl[t,h,x] = popt[0]
				else :
					print "++++ Erreur "
					exit()

######################
# Plot
######################
zeros=sp.ai_zeros(3)
xairy=np.linspace(-3,3,1000)
def p(h,chi,zeros):
        nominateur = np.power(sp.airy(abs(h)/chi+zeros[0])[0],2)
        denominateur = 2*integrate.quad(lambda x : np.power(sp.airy(abs(x)/chi+zeros[0])[0],2),0,np.inf)[0]
        return nominateur/denominateur
chi=np.power(integrate.quad(lambda x : np.power(sp.airy(abs(x)+zeros[1][0])[0],2),0,np.inf)[0]/integrate.quad(lambda x : np.power(x*sp.airy(abs(x)+zeros[1][0])[0],2),0,np.inf)[0],0.5)
valeur_airy = p(xairy,chi,zeros[1])

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
		valeur_moy=np.add(valeur[0:LX/2],valeur[LX/2:LX][::-1])
		valeur_moy=tuple([i/2 for i in valeur_moy])#-min(valeur_moy)
		valeur_moy -= min(valeur_moy)
		x=x[0:LX/2]
		moments = calcul_moment(x,valeur_moy)
                x=(x-moments[1])/moments[2];
                valeur_moy = np.array(valeur_moy)*moments[2]/moments[0]
		plt.plot(x,valeur_moy,label="t="+str(t)+" , h="+str(h))

	plt.plot(xairy,valeur_airy,".", label="airy", color="black")
	plt.plot(xairy,np.exp(-xairy*xairy/2)/np.power(2*np.pi,0.5),"-+",label="gaussienne",color="grey")
	plt.legend(prop={'size':10})
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
		valeur_moy=np.add(valeur[0:LX/2],valeur[LX/2:LX][::-1])
		valeur_moy=tuple([i/2 for i in valeur_moy])
	#	x,valeur = zip(*sorted(zip(x,valeur)))
		plt.plot(x[0:LX/2],valeur_moy,label="t="+str(t)+" , h="+str(h))
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
			popt,pcov=curve_fit(fit_tanh,data[:LX/2,0],data[:LX/2,1],p0=[3,15,1,0])	
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
	long_max=[]
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
		valeur_moy=np.add(valeur[0:LX/2],valeur[LX/2:LX][::-1])
		valeur_moy=tuple([i/2 for i in valeur_moy])
		popt,pcov=curve_fit(fit_gauss,x[0:LX/2],valeur_moy[0:LX/2],bounds=([0,10,0,-1],[8,20,8,1]))
		long_max=np.append(np.amax(valeur_moy[10:20]),long_max)
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
	plt.plot(champ,long_max,'+',label="t="+str(t)+" from maxs ",color=colors[t])
	plt.plot(ch_mag,l_mag,'--',label="t="+str(t)+" from mag profile",color=colors[t])
plt.legend()
plt.xlabel('Magnetic field')
plt.ylabel('X-Correlation length')
plt.title('X-correlation length after a gaussian fit of the Y-correlation length ')
plt.savefig(directory+'longueur_correl.pdf')
#plt.show()
