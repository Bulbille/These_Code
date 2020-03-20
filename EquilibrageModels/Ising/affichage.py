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

S=np.loadtxt(sys.argv[1])
M=int(S[-1][1])+1
N=int(S[-1][0])+1
s=np.zeros((M,N))

for k in range(len(S)):
    s[int(S[k,1]),int(S[k,0])]=int(S[k,2])

plt.figure()
plt.pcolormesh(s,rasterized=True,cmap="Pastel1")
plt.axis("image")
plt.axis("off")

#plt.savefig('t-'+str(temp)+'.pdf', bbox_inches='tight')
plt.show()
plt.close()

