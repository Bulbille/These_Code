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

J = 1.3 ; TC = 2/np.log(1.+np.power(2,0.5)) ; BETA = 1/TC;

directory=sys.argv[1]
taux = 0
if taux < 10 :
    regex_chain="X[\d]+_H1\.0e-0(\d)_F0\.0"+str(taux)
else :
    regex_chain="X[\d]+_H1\.0e-0(\d)_F0\."+str(taux)
colors=['red','blue','green','rosybrown','sandybrown','indigo','darkmagenta','gold','olivedrab'] ; nb_colors = len(colors)

###############################
# Mise en mémoire des données #
###############################
hspace = []
big_data = {}
for file in sorted(os.listdir(directory)):
    regex = re.search(regex_chain,file)
    if regex == None : 
            continue
    res = re.findall(regex_chain,file)
    if float(res[0]) not in hspace:
        hspace = np.append(hspace,int(res[0]))
        data=np.loadtxt(directory+file)
        big_data[float(res[0])] = data

hspace  = sorted(hspace)

############################
# Calcul de l'énergie libre
############################
def linear(x,a,b,p):
    return a*np.power(x,p)+b
def expo(x,a,b,p):
    return a*np.exp(p*x)+b

if taux == 0 :
    energie = plt.subplot(311)
    variance = plt.subplot(312)
    skew = plt.subplot(313)
else :
    drive = plt.subplot(221)
    energie = plt.subplot(222)
    variance = plt.subplot(223)
    skew = plt.subplot(224)
mm = 30

for nb,B in enumerate(hspace):
    print nb,B
    if nb != 1 : 
        continue
    if taux > 0 :
        drive.plot(big_data[B][:,0],big_data[B][:,1],color=colors[nb],    label='$B = 10^{-'+str(int(B))+'}$')
    variance.plot(big_data[B][:,0],big_data[B][:,2],color=colors[nb], label='$B = 10^{-'+str(int(B))+'}$')
    skew.plot(big_data[B][:,0],big_data[B][:,3],color=colors[nb],     label='$B = 10^{-'+str(int(B))+'}$')
    energie.plot(big_data[B][:,0],big_data[B][:,4],color=colors[nb],  label='$B = 10^{-'+str(int(B))+'}$')

    try :
        if taux > 0 :
            popt,pcov = curve_fit(linear,big_data[B][:mm,0],big_data[B][:mm,1],p0=[1,big_data[B][0,2],4])
            drive.plot(big_data[B][:mm,0],linear(big_data[B][:mm,0],*popt),'+',color=colors[nb])
            print "Mag", B,popt
        popt,pcov = curve_fit(linear,big_data[B][:mm,0],big_data[B][:mm,2],p0=[1,big_data[B][0,2],4])
        variance.plot(big_data[B][:mm,0],linear(big_data[B][:mm,0],*popt),'+',color=colors[nb])
        print "Var", B,popt
        popt,pcov = curve_fit(linear,big_data[B][:mm,0],big_data[B][:mm,4],p0=[1,big_data[B][0,2],4])
        energie.plot(big_data[B][:mm,0],linear(big_data[B][:mm,0],*popt),'+',color=colors[nb])
        print "Ene", B, popt
    except :
        pass

if taux > 0 :
    drive.set_xlabel('$f$')
    drive.set_ylabel('$\\bar{B}$')
    drive.legend(loc='upper left')
variance.set_xlabel('$f$')
energie.set_xlabel('$f$')
skew.set_xlabel('$f$')

variance.legend(loc='lower right')
energie.legend(loc='upper left')
skew.legend(loc='upper right')

variance.set_ylabel('$\sigma = <h^2>-<h>^2$')
energie.set_ylabel('$E=\sum |h_i-h_{i+1}|$')
skew.set_ylabel('$\gamma = <h-\\bar{h}>^3/\sigma^2$')

######ligne pointillé pour montrer 2J et 4J
#x1,y1 = [2*J,2*J],[0,5]
#x2,y2 = [4*J,4*J],y1
#energie.plot(x1,y1,'-.',x2,y2,'-.',color='black')
#variance.plot(x1,y1,'-.',x2,y2,'-.',color='black')
#skew.plot(x1,y1,'-.',x2,y2,'-.',color='black')

#variance.set_ylim([2.5,4.5])
#energie.set_ylim([1,8])
#skew.set_ylim([-6,-1])
variance.set_xlim([0,6])
energie.set_xlim([0,6])
skew.set_xlim([0,6])


plt.savefig(str(directory)+'/taux'+str(taux)+'.pdf')
plt.show()
plt.close()


