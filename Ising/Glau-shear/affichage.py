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
regex_chain="T-([\d.]+)(?:-\d)?(?:-t[\d.]*)?(?:-h-([\d.]+))?"
res = re.findall(regex_chain,sys.argv[1])

S=np.loadtxt(sys.argv[1])
M=int(S[-1][1])+1
N=int(S[-1][0])+1
s=np.zeros((M,N))

for k in range(len(S)):
		s[int(S[k,1]),int(S[k,0])]=int(S[k,2])

plt.figure()
plt.pcolormesh(s,rasterized=True,cmap="Pastel1")
plt.axis("image")
plt.xlabel("X")
plt.ylabel("Y")
print res
try:
    plt.title("Typical configuration at $T="+res[0][0]+"$ for $B=" + res[0][1]+ "$")
except :
    plt.title("Typical configuration at $T="+res[0][0]+"$ for $B=0$")

#try :
#	interface=np.loadtxt(sys.argv[1]+"-interface")
#	plt.plot(interface,'o')
#except :
#	pass

plt.savefig(sys.argv[1]+'-plot.png')
plt.show()
plt.close()

