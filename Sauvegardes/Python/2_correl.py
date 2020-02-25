#!/usr/bin/python
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

fontsize=15
plt.rc("font",size=fontsize)
plt.rc("font",family='serif')
plt.rc("text",usetex=True)

r0_255 = lambda: random.randint(0,255)
np.seterr(all='warn')

directory=sys.argv[1]
LX=40
LY=40
lx=np.arange(0,LX,1)
ly=np.arange(0,LY,1)
regex_chain="ttc_([\d.]+)(?:-\d)?(?:-t[\d.]*)?(?:-h-([\d.]+))?"

def calcul_moment(x,y):
        moment=[]
        x_bin=abs(x[1]-x[0])
	y=y-min(y)
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
colors=["red","blue","orange","green","black", "yellow", "purple", "pink", "salmon","red","blue","orange","green","black", "yellow", "purple", "pink", "salmon","red","blue","orange","green","black", "yellow", "purple", "pink", "salmon","red","blue","orange","green","black", "yellow", "purple", "pink", "salmon","red","blue","orange","green","black", "yellow", "purple", "pink", "salmon","red","blue","orange","green","black", "yellow", "purple", "pink", "salmon","red","blue","orange","green","black", "yellow", "purple", "pink", "salmon","red","blue","orange","green","black", "yellow", "purple", "pink", "salmon" ]
fmt=[":","--", ":","-"]

print "++++++++",temp,champs

############################
# Moyenne des données
############################
big_data = {}
nb = {}
big_mag = {}
big_energie = {}

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
		else :
			big_mag[t,h] = np.add(data,big_mag[t,h])
		
	elif re.search(r"energie",file) != None :
		res = re.findall(regex_chain,file)	
		t=float(res[0][0]); h=0
		if bool(res[0][1]):
			h=float(res[0][1])
		data=np.loadtxt(directory+file)
		if (t,h) not in big_energie:
			big_energie[t,h] = data
		else :
			big_energie[t,h] = np.add(data,big_energie[t,h])
	else : 
		continue



for points in big_data :
	big_data[points] /= nb[points]
	big_mag[points] /= nb[points]
	big_energie[points] /= nb[points]
####################################
# Profils de corrélation 
####################################
def fit_exp(x,l,s,a,h):
	#return -a*np.sin((x-s)/l)+h
	#return -a*np.tanh((x-s)/l)+h
	return a*np.exp(-(x-s)/l)+h

correl = {}

for t in temp:
	for h in champs :
		data=big_data[round(t,3),round(h,3)]
		for y in ly :
			datax=data[:][data[:,1]==y]
#			datax=datax[:][datax[:,1]<=LY/2]
#			print datax[:,0]
			popt=[]
#			print datax
#			plt.plot(datax[:,0],datax[:,1])
			try:
				popt,pcov=curve_fit(fit_exp,datax[:,0],datax[:,2]/np.amax(datax[:,2]),bounds=([0,-0.1,0,-1],[19.,0.1,1.1,1]))
			except :
				popt=np.append(popt,0)
			finally :
				if (t,h,y) not in correl:
					correl[t,h,y] = popt[0]
				else :
					print "++++ Erreur "
					exit()
######################
# Plot
######################
zeros=sp.ai_zeros(3)
xairy=np.linspace(-3,3,100)
def p(h,chi,zeros):
        nominateur = np.power(sp.airy(abs(h)/chi+zeros[0])[0],2)
        denominateur = 2*integrate.quad(lambda x : np.power(sp.airy(abs(x)/chi+zeros[0])[0],2),0,np.inf)[0]
        return nominateur/denominateur
chi=np.power(integrate.quad(lambda x : np.power(sp.airy(abs(x)+zeros[1][0])[0],2),0,np.inf)[0]/integrate.quad(lambda x : np.power(x*sp.airy(abs(x)+zeros[1][0])[0],2),0,np.inf)[0],0.5)
valeur_airy = p(xairy,chi,zeros[1])

for i,t in enumerate(temp) :
	for j,h in enumerate(champs) :
		label=True
		x=[]
		valeur=[]
		for points  in correl:
			if abs(float(points[0])-t)>0.0001 or abs(float(points[1])-h)>0.0001:
			#On n'est pas à la bonne température/champ mag	
				continue
			if points[2] > 7 and points[2] < LY - 7:
				x=np.append(x,points[2])
				valeur=np.append(valeur,correl[points])

		x,valeur = zip(*sorted(zip(x,valeur)))
		valeur = valeur-min(valeur)
		moments = calcul_moment(x,valeur)
		x=(x-moments[1])/moments[2]
		valeur = np.multiply(valeur,moments[2]/moments[0])
		if i == 0:
			plt.plot(x,valeur,fmt[i],label="$B=$"+str(h),color=colors[j])
		else :
			plt.plot(x,valeur,fmt[i],color=colors[j])

plt.plot(xairy,valeur_airy, label="airy", color="black")
plt.plot(xairy,np.exp(-xairy*xairy/2)/np.power(2*np.pi,0.5),label="gaussienne",color="grey")
plt.legend(loc=8,prop={'size':fontsize})
plt.xlabel('y')
plt.yscale('log')
plt.axis([-3.5,3.5,1e-3,5e-1])
plt.ylabel('$\\xi_\parallel$')
plt.title('')
plt.savefig(directory+'lc.pdf')
plt.show()
plt.close()
#for h in champs :
#	for t in temp :
#		x=[]
#		valeur=[]
#		for points  in correl:
#			if abs(float(points[0])-t)>0.0001 or abs(float(points[1])-h)>0.0001:
#			#On n'est pas à la bonne température/champ mag	
#				continue
#			x=np.append(x,points[2])
#			valeur=np.append(valeur,correl[points])
#		x,valeur = zip(*sorted(zip(x,valeur)))
#		moments = calcul_moment(x,valeur)
#                x=(x-moments[1])/moments[2]
#                valeur = np.multiply(valeur,moments[2]/moments[0])
#		plt.plot(x,valeur,label="t="+str(t)+" , h="+str(h))
#	plt.legend(loc='upper center',prop={'size':10})
#	plt.xlabel('X')
#	plt.ylabel('Y-Correlation length')
#	plt.title('Y-correlation length')
#	plt.savefig(directory+'lc_h'+str(h)+'.pdf')
##	plt.show()
#	plt.close()


####################################
# Fit du profil de magnétisation et moments de l'énergie
####################################
mag = {}
energie={}

def fit_tanh(x,l,s,a,h):
	return -a*np.tanh((x-s)/l)+h

for i,t in enumerate(temp) :
	for j,h in enumerate(champs) :
		datam=big_mag[t,h]
		datae=big_energie[t,h]
		popt=[]
		energie[t,h] = calcul_moment(datae[5:LY-5,0],datae[5:LY-5,1])[2]

		if i == 0:
			plt.figure(0)
			plt.plot(datam[:,0],datam[:,1],fmt[i],label="$B=$"+str(h), color=colors[j])
			plt.figure(1)
			plt.plot(datae[5:LY-4,0],datae[5:LY-4,1],fmt[i],label="$B=$"+str(h), color=colors[j])
		else :
			plt.figure(0)
			plt.plot(datam[:,0],datam[:,1],fmt[i],color=colors[j])
			plt.figure(1)
			plt.plot(datae[5:LY-4,0],datae[5:LY-4,1],fmt[i],color=colors[j])
		try:
			popt,pcov=curve_fit(fit_tanh,datam[:LX/2,0],datam[:LX/2,1],bounds=([0,LY/2-10,0.5,-1],[10,LY/2+10,1.5,1]))
		except RuntimeError :
			popt = np.append(popt,0)
		finally :
			if (t,h) not in mag:
				mag[t,h] = popt[0]
			else :
				mag[t,h] += popt[0]
plt.legend(prop={'size':fontsize})
plt.figure(0)
plt.legend(prop={'size':fontsize})
plt.savefig(directory+"mag.pdf")
plt.close()
plt.figure(1)
#plt.show()
plt.savefig(directory+"energie.pdf")
plt.close()

####################################
# Fit du profil en gaussienne
####################################
def fit_gauss(x,l,s,a,h):
	return a*np.exp(-(x-s)**2/(2*l**2))+h

for i,t in enumerate(temp) :
	champ=[]
	longueurs=[]
	long_max=[]
	for h in champs :
		label=True
		x=[]
		valeur=[]
		for points  in correl:
			if abs(float(points[0])-t)>0.0001 or abs(float(points[1])-h)>0.0001:
				continue
			if points[2] > 7 and points[2] < LY - 7:
				x=np.append(x,points[2])
				valeur=np.append(valeur,correl[points])

		x,valeur = zip(*sorted(zip(x,valeur)))
		valeur = valeur-min(valeur)
		moments = calcul_moment(x,valeur)

		champ=np.append(champ,h)
		longueurs=np.append(longueurs,np.abs(moments[2]))
	champ,longueurs = zip(*sorted(zip(champ,longueurs)))

	ch_mag=[] ; l_mag=[];
	for points in mag :
		if points[0] != t: 
			continue
		ch_mag=np.append(ch_mag,points[1])
		l_mag=np.append(l_mag,mag[points])
		ch_mag,l_mag = zip(*sorted(zip(ch_mag,l_mag)))
	ch_e=[] ; l_e = []
	for points in energie :
		if points[0] != t: 
			continue
		ch_e=np.append(ch_e,points[1])
		l_e=np.append(l_e,energie[points])
		ch_e,l_e = zip(*sorted(zip(ch_e,l_e)))
			
#	print "longueurs :",longueurs
#	print "mag : ", l_mag
#	print "energie : ",l_e
	plt.plot(champ,longueurs,'--',color=colors[i])
	plt.plot(ch_mag,l_mag,':',color=colors[i])
	plt.plot(ch_e,l_e,'-.',color=colors[i])
plt.axis([0.02,0.16, 1,10])
plt.xlabel('Magnetic field')
plt.ylabel('$\\xi_\\bot$')
plt.title('')
plt.savefig(directory+'longueur_correl.pdf')
plt.show()
