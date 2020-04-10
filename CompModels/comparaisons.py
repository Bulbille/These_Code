#!/usr/bin/python
# -*- coding: utf8


import numpy as np
import matplotlib.pyplot as plt
fontsize= 12
plt.rc("font",size=fontsize)
plt.rc("font",family='serif')
plt.rc("text",usetex=True) #Latex
import sys
import os
import re

#####################################
# Affiche l'état d'un système de 
# de taille N*M 
#####################################
LY  = 25.5
J   = 1
TC = 2/np.log(1+2**0.5)


def translation(x):
    print(2*J*(2-LY)/LY)
    return 2*J*(1-LY)/LY+2*x

dossier="bon"
try :
    ising =np.loadtxt('Ising/'+dossier+'/dataIsing')
    plt.plot(ising[:,0]/TC,ising[:,2],label="Ising")
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
try :
    pop =np.loadtxt('POP/'+dossier+'/dataPOP')
    plt.plot(pop[:,0],translation(pop[:,2]),label="POP")
except :
    pass
plt.legend()
plt.xlabel('$T$')
plt.ylabel('$\langle \mathcal{H} \\rangle$')

plt.savefig('comparaison-modeles.pdf')
plt.show()
plt.close()

