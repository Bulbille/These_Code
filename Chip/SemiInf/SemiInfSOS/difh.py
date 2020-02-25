#!/usr/bin/python
# -*- coding: utf8

#########################################
# Modèle d'interface en histogramme
#########################################
import matplotlib.pyplot as plt
fontsize=9
labelsize=16
plt.rc("font",size=fontsize)
plt.rc("font",family='serif')

import scipy as sy
from scipy.optimize import curve_fit
import numpy as np
from numpy import linalg as LA
import os #pour les chemins de fichiers
import sys #pour l'exit et l'argument
import re

TC = 2/np.log(1.+np.power(2,0.5));BETA = 1/TC
J = 1

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

if taux == 0 :
    energie = plt.subplot(311)
    variance = plt.subplot(312)
    skew = plt.subplot(313)
else :
    drive = plt.subplot(221)
    energie = plt.subplot(222)
    variance = plt.subplot(223)
    skew = plt.subplot(224)
derivene = energie.twinx()
mm = 10

for nb,B in enumerate(hspace):
    print "**************" ,nb,B
    fspace = big_data[B][:,0]
    donne  = np.where(fspace<6)
    fspace = fspace[donne]
    totmag = big_data[B][donne,1][0]
    totvar = big_data[B][donne,2][0]
    totske = big_data[B][donne,3][0]
    totene = big_data[B][donne,4][0]
#    if nb >= 1 : 
#        continue
    if taux > 0 :
        drive.plot(fspace,totmag,color=colors[nb],    label='$B = 10^{-'+str(int(B))+'}$')
    variance.plot(fspace,totvar**2/totmag**2,color=colors[nb], label='$B = 10^{-'+str(int(B))+'}$')
    skew.plot(fspace,totske**3/totmag**3,color=colors[nb],     label='$B = 10^{-'+str(int(B))+'}$')
    energie.plot(fspace,totene/totene[0],color=colors[nb],  label='$B = 10^{-'+str(int(B))+'}$')
    derivene.plot(fspace,np.gradient(totene/totene[0]),'-.',color=colors[nb],  label='$B = 10^{-'+str(int(B))+'}$')
    try :
        if taux > 0 :
            popt,pcov = curve_fit(linear,fspace[:mm],totmag[:mm],p0=[1,1,2])
            drive.plot(fspace[:mm],linear(fspace[:mm],*popt),'+',color=colors[nb])
            print "Mag", B,popt
        popt,pcov = curve_fit(linear,fspace[:mm],totvar[:mm]**2/totmag[:mm]**2,p0=[1,1,2])
        variance.plot(fspace[:mm],linear(fspace[:mm],*popt),'+',color=colors[nb])
        print "Var", B,popt
        popt,pcov = curve_fit(linear,fspace[:mm],totene[:mm]/totene[0],p0=[1,1,2])
        energie.plot(fspace[:mm],linear[:mm](fspace,*popt),'+',color=colors[nb])
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

variance.legend(loc='upper left')
energie.legend(loc='upper left')
skew.legend(loc='lower right')

energie.set_ylabel('$E/E_0$',fontsize=labelsize)
variance.set_ylabel('$\sigma^2/H^2$',fontsize=labelsize)
skew.set_ylabel('$\gamma^3/H^3$',fontsize=labelsize)

#ligne pointillé pour montrer 2J et 4J
energie.axvline(2*J,linestyle='-.',color='black')
energie.axvline(4*J,linestyle='-.',color='black')
variance.axvline(2*J,linestyle='-.',color='black')
variance.axvline(4*J,linestyle='-.',color='black')
skew.axvline(2*J,linestyle='-.',color='black')
skew.axvline(4*J,linestyle='-.',color='black')
#variance.plot(x1,y1,'-.',x2,y2,'-.',color='black')
#skew.plot(x1,y1,'-.',x2,y2,'-.',color='black')

#variance.set_ylim([2.5,4.5])
#energie.set_ylim([1,8])
#skew.set_ylim([-6,-1])


plt.savefig(str(directory)+'/taux'+str(taux)+'.pdf')
plt.show()
plt.close()


