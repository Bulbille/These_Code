#!/usr/bin/python
# -*- coding: utf8


import numpy as np
import matplotlib.pyplot as plt
fontsize= 12
plt.rc("font",size=fontsize)
plt.rc("font",family='serif')
#plt.rc("text",usetex=True) #Latex
import sys
import os
import re

#####################################
# Affiche l'état d'un système de 
# de taille N*M 
#####################################
LY  = 24
J   = 1


def translation(x):
    print(2*J*(2-LY)/LY)
    return 2*J*(2-LY)/LY+2*x

dossier="bon"
try :
    ising =np.loadtxt('Ising/'+dossier+'/dataIsing')
    plt.plot(ising[:,0],ising[:,2],label="Ising")
except :
    pass
try :
    sos =np.loadtxt('SOS/'+dossier+'/dataSOS')
    plt.plot(sos[:,0],translation(sos[:,3]),label="SOS")
except :
    pass
try :
    rsos =np.loadtxt('RSOS/'+dossier+'/dataRSOS')
    plt.plot(rsos[:,0],translation(rsos[:,3]),label="RSOS")
except :
    pass
plt.legend()
#plt.xlabel('$\\frac{T}{T_{C,2D}}$')
#plt.ylabel('Énergie du système')


plt.show()
plt.close()

